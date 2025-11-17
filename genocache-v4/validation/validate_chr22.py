#!/usr/bin/env python3
"""
Validate trained model on chr22 test reads

This script:
1. Loads trained chr22 model
2. Encodes test reads (512bp)
3. Encodes chr22 genome with stride=32 (CRITICAL!)
4. Builds FAISS IVFPQ index
5. Searches and measures accuracy

Target: >80% accuracy at ±1kb tolerance
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
import faiss
import json
from pathlib import Path
from tqdm import tqdm
import time

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder


def load_genome(fasta_path: Path, chr_name: str = "NC_000022.11") -> str:
    """Load specific chromosome from genome"""
    print(f"Loading {chr_name} from {fasta_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr is not None:
                    genome[current_chr] = ''.join(current_seq).upper()
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_chr is not None:
            genome[current_chr] = ''.join(current_seq).upper()
    
    # Try different naming conventions
    for name in [chr_name, 'chr22', '22', 'NC_000022.11']:
        if name in genome:
            print(f"✅ Found as '{name}': {len(genome[name]):,} bp")
            return genome[name]
    
    raise ValueError(f"Could not find {chr_name}!")


def encode_genome_with_stride(
    model: GenoCacheEncoder,
    genome_seq: str,
    window_size: int = 512,
    stride: int = 32,
    device: str = 'cuda',
    batch_size: int = 256
) -> tuple:
    """
    Encode genome with sparse indexing (stride=32)
    
    CRITICAL: This enables 32× memory reduction via translation continuity
    
    Args:
        model: Trained encoder
        genome_seq: Chromosome sequence
        window_size: Window size (512bp to match training)
        stride: Stride for sparse indexing (32bp)
        device: Device
        batch_size: Batch size for encoding
    
    Returns:
        (embeddings, positions): Embeddings array and their positions
    """
    model.eval()
    
    print(f"Encoding genome with stride={stride}...")
    print(f"  Genome length: {len(genome_seq):,} bp")
    print(f"  Window size: {window_size} bp")
    print(f"  Stride: {stride} bp")
    
    # Calculate number of windows
    num_windows = (len(genome_seq) - window_size) // stride + 1
    print(f"  Number of windows: {num_windows:,}")
    print(f"  Compared to dense: {len(genome_seq) - window_size + 1:,} (saved {stride}×)")
    
    embeddings_list = []
    positions_list = []
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        """Convert DNA sequence to index tensor"""
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    batch_seqs = []
    batch_positions = []
    
    with torch.no_grad():
        for pos in tqdm(range(0, len(genome_seq) - window_size + 1, stride), desc="Encoding"):
            window = genome_seq[pos:pos + window_size]
            
            # Skip high-N regions
            if window.count('N') > window_size * 0.1:
                continue
            
            batch_seqs.append(seq_to_tensor(window))
            batch_positions.append(pos)
            
            # Process batch
            if len(batch_seqs) >= batch_size:
                batch_tensor = torch.stack(batch_seqs).to(device)
                embeddings = model(batch_tensor)
                embeddings = F.normalize(embeddings, p=2, dim=1)
                
                embeddings_list.append(embeddings.cpu().numpy())
                positions_list.extend(batch_positions)
                
                batch_seqs = []
                batch_positions = []
        
        # Process remaining
        if batch_seqs:
            batch_tensor = torch.stack(batch_seqs).to(device)
            embeddings = model(batch_tensor)
            embeddings = F.normalize(embeddings, p=2, dim=1)
            
            embeddings_list.append(embeddings.cpu().numpy())
            positions_list.extend(batch_positions)
    
    # Concatenate all
    all_embeddings = np.vstack(embeddings_list)
    all_positions = np.array(positions_list)
    
    print(f"✅ Encoded {len(all_positions):,} windows")
    print(f"  Embedding shape: {all_embeddings.shape}")
    print(f"  Memory: {all_embeddings.nbytes / 1e6:.1f} MB")
    
    return all_embeddings, all_positions


def build_faiss_index(embeddings: np.ndarray, nprobe: int = 16) -> faiss.IndexIVFPQ:
    """
    Build FAISS IVFPQ index
    
    Args:
        embeddings: [N, D] embeddings
        nprobe: Number of clusters to probe (higher = better accuracy, slower)
    
    Returns:
        FAISS index
    """
    N, D = embeddings.shape
    print(f"\nBuilding FAISS index...")
    print(f"  Embeddings: {N:,} × {D}D")
    
    # Configuration following NeuralAligner paper
    nlist = int(np.sqrt(N))  # Number of clusters
    nlist = max(nlist, 100)  # Minimum 100 clusters
    nlist = min(nlist, 4096)  # Maximum 4096 clusters
    
    print(f"  Index type: IVFPQ")
    print(f"  Number of clusters (nlist): {nlist}")
    print(f"  nprobe: {nprobe}")
    
    # Create index
    quantizer = faiss.IndexFlatIP(D)  # Inner product (for cosine similarity)
    index = faiss.IndexIVFPQ(quantizer, D, nlist, 16, 8)
    
    # Train and add
    print("  Training index...")
    index.train(embeddings.astype(np.float32))
    
    print("  Adding embeddings...")
    index.add(embeddings.astype(np.float32))
    
    # Set search parameters
    index.nprobe = nprobe
    
    print(f"✅ Index built: {index.ntotal:,} embeddings")
    
    return index


def load_test_reads(data_dir: Path) -> tuple:
    """Load test reads and ground truth"""
    print("Loading test reads...")
    
    # Load sequences
    seq_file = data_dir / "test_reads.fasta"
    reads = []
    with open(seq_file, 'r') as f:
        seq = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq is not None:
                    reads.append(seq)
                seq = None
            else:
                seq = line
        if seq is not None:
            reads.append(seq)
    
    # Load ground truth
    truth_file = data_dir / "ground_truth.json"
    with open(truth_file, 'r') as f:
        ground_truth = json.load(f)
    
    print(f"✅ Loaded {len(reads)} test reads")
    
    return reads, ground_truth


def encode_reads(model: GenoCacheEncoder, reads: list, device: str = 'cuda', batch_size: int = 64):
    """Encode test reads"""
    model.eval()
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, len(reads), batch_size), desc="Encoding reads"):
            batch = reads[i:i+batch_size]
            batch_tensors = torch.stack([seq_to_tensor(seq) for seq in batch]).to(device)
            
            batch_emb = model(batch_tensors)
            batch_emb = F.normalize(batch_emb, p=2, dim=1)
            
            embeddings.append(batch_emb.cpu().numpy())
    
    return np.vstack(embeddings)


def validate(
    index: faiss.IndexIVFPQ,
    positions: np.ndarray,
    read_embeddings: np.ndarray,
    ground_truth: list,
    k: int = 32,
    tolerance: int = 1000
) -> dict:
    """
    Run validation and measure accuracy
    
    Args:
        index: FAISS index
        positions: Genome positions for each embedding
        read_embeddings: Test read embeddings
        ground_truth: Ground truth positions
        k: Number of top matches to retrieve
        tolerance: Tolerance in bp for correct match (±1000bp)
    
    Returns:
        metrics: Accuracy metrics
    """
    print(f"\nValidating...")
    print(f"  Reads: {len(read_embeddings)}")
    print(f"  Top-k: {k}")
    print(f"  Tolerance: ±{tolerance} bp")
    
    correct_1kb = 0
    correct_100bp = 0
    correct_exact = 0
    
    top1_errors = []
    best_k_errors = []
    
    for i in tqdm(range(len(read_embeddings)), desc="Searching"):
        query_emb = read_embeddings[i:i+1].astype(np.float32)
        true_pos = ground_truth[i]['true_start']
        
        # Search
        distances, indices = index.search(query_emb, k)
        
        # Get predicted positions
        pred_positions = positions[indices[0]]
        
        # Check top-1 accuracy
        top1_pos = pred_positions[0]
        top1_error = abs(top1_pos - true_pos)
        top1_errors.append(top1_error)
        
        # Check best-of-k accuracy
        best_error = min(abs(p - true_pos) for p in pred_positions)
        best_k_errors.append(best_error)
        
        # Count correct
        if best_error <= 1000:
            correct_1kb += 1
        if best_error <= 100:
            correct_100bp += 1
        if best_error == 0:
            correct_exact += 1
    
    # Calculate metrics
    n_reads = len(read_embeddings)
    metrics = {
        'total_reads': n_reads,
        'accuracy_1kb': 100 * correct_1kb / n_reads,
        'accuracy_100bp': 100 * correct_100bp / n_reads,
        'accuracy_exact': 100 * correct_exact / n_reads,
        'top1_median_error': np.median(top1_errors),
        'top1_mean_error': np.mean(top1_errors),
        'bestk_median_error': np.median(best_k_errors),
        'bestk_mean_error': np.mean(best_k_errors),
    }
    
    return metrics


def main():
    print("=" * 80)
    print("GenoCache V4 - chr22 Validation")
    print("=" * 80)
    print()
    
    # Configuration
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    CHECKPOINT_DIR = Path("/home/nebius/genocache/genocache-v4/models/checkpoints")
    DATA_DIR = Path("/home/nebius/genocache/genocache-v4/validation/data")
    OUTPUT_DIR = Path("/home/nebius/genocache/genocache-v4/validation/results")
    
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    # Find best checkpoint
    checkpoints = list(CHECKPOINT_DIR.glob("chr22_best_*.pt"))
    if not checkpoints:
        # Try any chr22 checkpoint
        checkpoints = list(CHECKPOINT_DIR.glob("chr22_epoch*.pt"))
    
    if not checkpoints:
        print("❌ No chr22 checkpoint found!")
        return
    
    # Use the latest one
    checkpoint_path = sorted(checkpoints)[-1]
    print(f"Checkpoint: {checkpoint_path.name}\n")
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}\n")
    
    # Load model
    print("Loading model...")
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    num_params = sum(p.numel() for p in model.parameters())
    print(f"✅ Model loaded: {num_params:,} parameters")
    
    # Load chr22
    chr22_seq = load_genome(GENOME_PATH)
    
    # Encode genome with stride=32 (CRITICAL!)
    start_time = time.time()
    genome_embeddings, genome_positions = encode_genome_with_stride(
        model, chr22_seq,
        window_size=512,
        stride=32,  # Sparse indexing!
        device=device,
        batch_size=256
    )
    genome_time = time.time() - start_time
    print(f"  Encoding time: {genome_time:.1f}s")
    
    # Build FAISS index
    start_time = time.time()
    index = build_faiss_index(genome_embeddings, nprobe=16)
    index_time = time.time() - start_time
    print(f"  Index building time: {index_time:.1f}s")
    
    # Load test reads
    reads, ground_truth = load_test_reads(DATA_DIR)
    
    # Encode test reads
    print("\nEncoding test reads...")
    start_time = time.time()
    read_embeddings = encode_reads(model, reads, device=device, batch_size=64)
    read_time = time.time() - start_time
    print(f"✅ Encoded {len(reads)} reads in {read_time:.1f}s ({read_time/len(reads)*1000:.1f}ms per read)")
    
    # Validate
    metrics = validate(
        index, genome_positions, read_embeddings, ground_truth,
        k=32, tolerance=1000
    )
    
    # Print results
    print("\n" + "=" * 80)
    print("VALIDATION RESULTS")
    print("=" * 80)
    print(f"Total reads: {metrics['total_reads']}")
    print()
    print(f"Accuracy @ ±1kb:   {metrics['accuracy_1kb']:.2f}%")
    print(f"Accuracy @ ±100bp: {metrics['accuracy_100bp']:.2f}%")
    print(f"Accuracy @ exact:  {metrics['accuracy_exact']:.2f}%")
    print()
    print(f"Top-1 median error: {metrics['top1_median_error']:.0f} bp")
    print(f"Top-1 mean error:   {metrics['top1_mean_error']:.0f} bp")
    print(f"Best-k median error: {metrics['bestk_median_error']:.0f} bp")
    print(f"Best-k mean error:   {metrics['bestk_mean_error']:.0f} bp")
    print()
    
    # Check target
    print("Target Check:")
    if metrics['accuracy_1kb'] >= 80:
        print(f"  ✅ PASSED! Accuracy {metrics['accuracy_1kb']:.2f}% >= 80%")
        print(f"  → Ready to scale to full genome!")
    elif metrics['accuracy_1kb'] >= 60:
        print(f"  ⚠️  MARGINAL: Accuracy {metrics['accuracy_1kb']:.2f}% (60-80%)")
        print(f"  → Analyze failures, may need curriculum learning")
    else:
        print(f"  ❌ FAILED: Accuracy {metrics['accuracy_1kb']:.2f}% < 60%")
        print(f"  → Debug embeddings, check model convergence")
    
    # Save results
    results_file = OUTPUT_DIR / f"chr22_validation_{checkpoint_path.stem}.json"
    with open(results_file, 'w') as f:
        json.dump({
            'checkpoint': str(checkpoint_path),
            'metrics': metrics,
            'timing': {
                'genome_encoding': genome_time,
                'index_building': index_time,
                'read_encoding': read_time
            }
        }, f, indent=2)
    
    print(f"\n✅ Results saved: {results_file}")


if __name__ == "__main__":
    main()
