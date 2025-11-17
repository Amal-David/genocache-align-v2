#!/usr/bin/env python3
"""
NAL-Aligned Index Builder - Exact NeuralAligner Protocol

Following NAL paper Section 3.2 and A.3:
- FAISS with IVFPQ (Inverted File + Product Quantization)
- Stride: 32 (NAL uses 16-32, we use 32 for 512bp seeds)
- nlist = sqrt(num_vectors)
- nprobe = 8-32 (we use 8 for speed, 32 for accuracy)
- K = 32 (top-K neighbors per seed)
- Distance: inner product
- Compression: PQ16×8 (128D → 24 bytes)
- Use encoder output DIRECTLY (no projection head!)

Reference: NeuralAligner paper Section A.3
"""

import sys
import torch
import faiss
import numpy as np
from pathlib import Path
from tqdm import tqdm
import argparse
from datetime import datetime

# Add paths
sys.path.insert(0, str(Path(__file__).parent))
from encoder_nal import NALEncoder


def load_genome(fasta_path: Path) -> dict:
    """Load genome from FASTA"""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous chromosome
                if current_chr is not None:
                    seq = ''.join(current_seq).upper()
                    # Filter out high-N regions (NAL does this)
                    if seq.count('N') / len(seq) < 0.1:
                        genome[current_chr] = seq
                
                # Start new chromosome
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last chromosome
        if current_chr is not None:
            seq = ''.join(current_seq).upper()
            if seq.count('N') / len(seq) < 0.1:
                genome[current_chr] = seq
    
    total_len = sum(len(s) for s in genome.values())
    print(f"  Loaded {len(genome)} chromosomes, {total_len:,} bp total")
    
    return genome


def extract_seeds(genome_dict, seed_len=512, stride=32):
    """
    Extract seeds with sliding window (NAL protocol)
    
    NAL: stride ≤ seed_len/8 (for 512bp, max stride is 64)
    We use stride=32 for good coverage
    
    Args:
        genome_dict: Dictionary of chromosome sequences
        seed_len: Seed length (NAL uses 256-512)
        stride: Step size (NAL uses 16-32)
    
    Returns:
        seeds: List of DNA sequences
        positions: List of (chr, pos) tuples
    """
    print(f"\nExtracting seeds (length={seed_len}, stride={stride})...")
    
    seeds = []
    positions = []
    
    for chr_name, seq in tqdm(genome_dict.items(), desc="Chromosomes"):
        # Slide window across chromosome
        for i in range(0, len(seq) - seed_len + 1, stride):
            seed = seq[i:i + seed_len]
            
            # Skip if too many N's (NAL does this)
            if seed.count('N') > seed_len * 0.1:
                continue
            
            seeds.append(seed)
            positions.append((chr_name, i))
    
    print(f"  Extracted {len(seeds):,} seeds")
    print(f"  Coverage: ~{len(seeds) * stride / sum(len(s) for s in genome_dict.values()) * 100:.1f}%")
    
    return seeds, positions


def encode_seeds(model, seeds, device='cuda', batch_size=1024, use_projection=False):
    """
    Encode seeds to embeddings (NAL: use encoder output directly!)
    
    CRITICAL: use_projection=False during indexing (NAL Section A.3)
    Projection head is ONLY for training, NOT for inference/indexing!
    
    Args:
        model: NAL encoder
        seeds: List of DNA sequences
        device: Device to use
        batch_size: Batch size for encoding
        use_projection: Whether to use projection head (should be False!)
    
    Returns:
        embeddings: (N, 128) array of embeddings
    """
    print(f"\nEncoding seeds (batch_size={batch_size}, projection={use_projection})...")
    
    model.eval()
    model.to(device)
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    embeddings = []
    
    # Process in batches
    num_batches = (len(seeds) + batch_size - 1) // batch_size
    
    with torch.no_grad():
        for i in tqdm(range(0, len(seeds), batch_size), total=num_batches, desc="Encoding"):
            batch_seeds = seeds[i:i + batch_size]
            
            # Convert to tensors
            batch_tensors = []
            for seed in batch_seeds:
                indices = [base_to_idx.get(b, 4) for b in seed]
                batch_tensors.append(torch.tensor(indices, dtype=torch.long))
            
            # Stack and move to device
            x = torch.stack(batch_tensors).to(device)
            
            # Encode (no projection head!)
            emb = model(x, use_projection=use_projection)
            
            # Move to CPU and convert to numpy
            embeddings.append(emb.cpu().numpy())
    
    embeddings = np.vstack(embeddings)
    print(f"  Encoded {embeddings.shape[0]:,} seeds → {embeddings.shape[1]}D")
    print(f"  Memory: {embeddings.nbytes / 1e9:.2f} GB (FP32)")
    
    return embeddings


def build_faiss_index(embeddings, nlist=None, nprobe=8, use_gpu=True):
    """
    Build FAISS index (NAL: IVFPQ with PQ16×8)
    
    NAL configuration (Section A.3):
    - Index type: IVFPQ (Inverted File + Product Quantization)
    - nlist = sqrt(num_vectors)
    - nprobe = 8-32
    - PQ: 16 subquantizers × 8 bits = 24 bytes per vector
    - Distance: inner product
    
    Args:
        embeddings: (N, D) array
        nlist: Number of clusters (default: sqrt(N))
        nprobe: Number of clusters to probe (8-32)
        use_gpu: Whether to use GPU
    
    Returns:
        index: FAISS index
    """
    print(f"\nBuilding FAISS index (IVFPQ)...")
    
    n_vectors, dim = embeddings.shape
    
    # NAL: nlist = sqrt(num_vectors)
    if nlist is None:
        nlist = int(np.sqrt(n_vectors))
    
    print(f"  Vectors: {n_vectors:,}")
    print(f"  Dimension: {dim}")
    print(f"  nlist (clusters): {nlist}")
    print(f"  nprobe: {nprobe}")
    
    # Create quantizer (NAL uses inner product)
    quantizer = faiss.IndexFlatIP(dim)
    
    # Create IVFPQ index
    # PQ16x8: 16 subquantizers, 8 bits each
    # 128D / 16 = 8D per subquantizer
    # Compression: 128 × 4 bytes (FP32) → 16 bytes (PQ) = 8× compression
    index = faiss.IndexIVFPQ(quantizer, dim, nlist, 16, 8)
    
    # Set distance metric to inner product (NAL uses this)
    index.metric_type = faiss.METRIC_INNER_PRODUCT
    
    if use_gpu:
        try:
            print("  Attempting to use GPU for training...")
            res = faiss.StandardGpuResources()
            index = faiss.index_cpu_to_gpu(res, 0, index)
            print("  GPU enabled successfully!")
        except (AttributeError, RuntimeError) as e:
            print(f"  GPU not available ({e}), using CPU instead...")
            use_gpu = False
    
    # Train index (required for IVFPQ)
    print("  Training index...")
    # Use subset for training if dataset is large
    train_size = min(n_vectors, 256 * nlist)
    train_indices = np.random.choice(n_vectors, train_size, replace=False)
    train_data = embeddings[train_indices].astype('float32')
    
    # Normalize for inner product (NAL normalizes embeddings)
    faiss.normalize_L2(train_data)
    
    index.train(train_data)
    
    # Add all vectors
    print("  Adding vectors...")
    # Normalize all embeddings
    embeddings_norm = embeddings.astype('float32')
    faiss.normalize_L2(embeddings_norm)
    
    # Add in batches to avoid memory issues
    batch_size = 100000
    for i in tqdm(range(0, n_vectors, batch_size), desc="Adding"):
        batch = embeddings_norm[i:i + batch_size]
        index.add(batch)
    
    # Set search parameters
    index.nprobe = nprobe
    
    # Move back to CPU for saving (if on GPU)
    if use_gpu:
        try:
            print("  Moving index to CPU for saving...")
            index = faiss.index_gpu_to_cpu(index)
        except:
            pass  # Already on CPU
    
    print(f"  Index built: {index.ntotal:,} vectors")
    
    # Calculate memory usage
    # PQ16x8: 16 bytes per vector + overhead
    memory_mb = index.ntotal * 16 / 1e6
    print(f"  Estimated memory: {memory_mb:.1f} MB (compressed)")
    
    return index


def save_index(index, positions, output_path, metadata=None):
    """
    Save FAISS index and position mapping
    
    Args:
        index: FAISS index
        positions: List of (chr, pos) tuples
        output_path: Output path for index
        metadata: Optional metadata dict
    """
    print(f"\nSaving index to {output_path}...")
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save FAISS index
    faiss.write_index(index, str(output_path))
    print(f"  Saved FAISS index: {output_path}")
    
    # Save position mapping
    pos_file = output_path.with_suffix('.positions.npz')
    chr_names = [p[0] for p in positions]
    chr_positions = np.array([p[1] for p in positions], dtype=np.int64)
    
    np.savez_compressed(
        pos_file,
        chr_names=chr_names,
        positions=chr_positions,
        metadata=metadata or {}
    )
    print(f"  Saved positions: {pos_file}")
    
    # Print summary
    print(f"\n✅ Index building complete!")
    print(f"  Total vectors: {index.ntotal:,}")
    print(f"  Index size: {output_path.stat().st_size / 1e6:.1f} MB")
    print(f"  Positions size: {pos_file.stat().st_size / 1e6:.1f} MB")


def main():
    parser = argparse.ArgumentParser(description='Build NAL-aligned FAISS index')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--model', required=True, help='Trained NAL model')
    parser.add_argument('--output', required=True, help='Output index file')
    parser.add_argument('--seed-len', type=int, default=512, help='Seed length (NAL: 256-512)')
    parser.add_argument('--stride', type=int, default=32, help='Stride (NAL: ≤seed_len/8)')
    parser.add_argument('--nprobe', type=int, default=8, help='Number of clusters to probe (8-32)')
    parser.add_argument('--batch-size', type=int, default=1024, help='Encoding batch size')
    parser.add_argument('--device', default='cuda', choices=['cuda', 'cpu'])
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("NAL-Aligned Index Builder")
    print("=" * 80)
    print(f"\nConfiguration (following NAL paper):")
    print(f"  Genome: {args.genome}")
    print(f"  Model: {args.model}")
    print(f"  Seed length: {args.seed_len} bp")
    print(f"  Stride: {args.stride} bp (NAL: ≤{args.seed_len//8})")
    print(f"  nprobe: {args.nprobe} (NAL: 8-32)")
    print(f"  Device: {args.device}")
    print()
    
    # Validate stride
    if args.stride > args.seed_len // 8:
        print(f"⚠️  WARNING: stride ({args.stride}) > seed_len/8 ({args.seed_len//8})")
        print(f"    NAL recommends stride ≤ seed_len/8 for good accuracy")
    
    # Load model
    print("Loading NAL model...")
    checkpoint = torch.load(args.model, map_location='cpu')
    model = NALEncoder(
        emb_dim=128,
        seed_len=args.seed_len,
        vocab_size=5,
        hidden_dim=128,
        num_layers=4
    )
    model.load_state_dict(checkpoint['model_state_dict'])
    print(f"  Loaded from epoch {checkpoint['epoch']}")
    print(f"  Val loss: {checkpoint['val_metrics']['loss']:.4f}")
    
    # Load genome
    genome_dict = load_genome(Path(args.genome))
    
    # Extract seeds
    seeds, positions = extract_seeds(
        genome_dict,
        seed_len=args.seed_len,
        stride=args.stride
    )
    
    # Encode seeds (NO projection head!)
    embeddings = encode_seeds(
        model,
        seeds,
        device=args.device,
        batch_size=args.batch_size,
        use_projection=False  # CRITICAL: NAL uses encoder output directly!
    )
    
    # Build index
    index = build_faiss_index(
        embeddings,
        nlist=None,  # Will be set to sqrt(n_vectors)
        nprobe=args.nprobe,
        use_gpu=(args.device == 'cuda')
    )
    
    # Save index
    metadata = {
        'seed_len': args.seed_len,
        'stride': args.stride,
        'nprobe': args.nprobe,
        'model_epoch': checkpoint['epoch'],
        'val_loss': checkpoint['val_metrics']['loss'],
        'build_date': datetime.now().isoformat(),
        'num_vectors': len(seeds),
        'genome_size': sum(len(s) for s in genome_dict.values())
    }
    
    save_index(index, positions, args.output, metadata)
    
    print("\n" + "=" * 80)
    print("Index ready for alignment!")
    print("=" * 80)


if __name__ == '__main__':
    main()
