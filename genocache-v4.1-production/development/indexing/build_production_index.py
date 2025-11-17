#!/usr/bin/env python3
"""
Build Production-Ready Full Genome Index

This builds the complete GenoCache index for production use:
1. Encode all primary chromosomes (chr1-22, X, Y)
2. Use stride=32 for sparse indexing
3. Build optimized FAISS IVFPQ index
4. Save index to disk for reuse
5. Benchmark index building time and memory
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
import faiss
import pickle
import time
from pathlib import Path
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder


def load_genome(fasta_path: Path, chromosomes: list = None) -> dict:
    """Load specified chromosomes from genome"""
    print(f"Loading genome from {fasta_path}...")
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
    
    # Filter to requested chromosomes if specified
    if chromosomes:
        genome = {k: v for k, v in genome.items() if k in chromosomes}
    
    return genome


def encode_genome(
    model: torch.nn.Module,
    genome: dict,
    window_size: int = 512,
    stride: int = 32,
    device: str = 'cuda',
    batch_size: int = 512
) -> tuple:
    """
    Encode genome with sparse indexing
    
    Returns:
        (embeddings, positions, chr_names): Arrays for index
    """
    model.eval()
    
    print(f"\nEncoding genome...")
    print(f"  Chromosomes: {len(genome)}")
    print(f"  Window size: {window_size} bp")
    print(f"  Stride: {stride} bp")
    
    # Calculate total windows
    total_windows = sum((len(seq) - window_size) // stride + 1 for seq in genome.values())
    total_bp = sum(len(seq) for seq in genome.values())
    print(f"  Total genome size: {total_bp:,} bp")
    print(f"  Total windows: {total_windows:,}")
    print(f"  Memory reduction: {stride}× vs dense indexing")
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    all_embeddings = []
    all_positions = []
    all_chr_names = []
    
    start_time = time.time()
    
    with torch.no_grad():
        for chr_name in tqdm(sorted(genome.keys()), desc="Encoding chromosomes"):
            chr_seq = genome[chr_name]
            chr_start_time = time.time()
            
            batch_seqs = []
            batch_positions = []
            chr_embeddings = 0
            
            for pos in range(0, len(chr_seq) - window_size + 1, stride):
                window = chr_seq[pos:pos + window_size]
                
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
                    
                    all_embeddings.append(embeddings.cpu().numpy())
                    all_positions.extend(batch_positions)
                    all_chr_names.extend([chr_name] * len(batch_positions))
                    
                    chr_embeddings += len(batch_positions)
                    batch_seqs = []
                    batch_positions = []
            
            # Process remaining
            if batch_seqs:
                batch_tensor = torch.stack(batch_seqs).to(device)
                embeddings = model(batch_tensor)
                embeddings = F.normalize(embeddings, p=2, dim=1)
                
                all_embeddings.append(embeddings.cpu().numpy())
                all_positions.extend(batch_positions)
                all_chr_names.extend([chr_name] * len(batch_positions))
                chr_embeddings += len(batch_positions)
            
            chr_time = time.time() - chr_start_time
            print(f"    {chr_name}: {len(chr_seq):,} bp → {chr_embeddings:,} embeddings ({chr_time:.1f}s)")
    
    # Concatenate all
    all_embeddings = np.vstack(all_embeddings)
    all_positions = np.array(all_positions, dtype=np.int64)
    all_chr_names = np.array(all_chr_names)
    
    encoding_time = time.time() - start_time
    
    print(f"\n✅ Genome encoding complete!")
    print(f"  Total embeddings: {len(all_positions):,}")
    print(f"  Embedding shape: {all_embeddings.shape}")
    print(f"  Memory: {all_embeddings.nbytes / 1e9:.2f} GB")
    print(f"  Encoding time: {encoding_time/60:.1f} minutes")
    print(f"  Speed: {total_bp / encoding_time / 1e6:.2f} Mbp/s")
    
    return all_embeddings, all_positions, all_chr_names


def build_faiss_index(embeddings: np.ndarray, nprobe: int = 64) -> faiss.Index:
    """Build optimized FAISS IVFPQ index for production"""
    N, D = embeddings.shape
    print(f"\nBuilding FAISS index...")
    print(f"  Embeddings: {N:,} × {D}D")
    
    # Optimize number of clusters
    nlist = min(int(np.sqrt(N)), 8192)  # Cap at 8192 for efficiency
    print(f"  Index type: IVFPQ (Inverted File + Product Quantization)")
    print(f"  Number of clusters (nlist): {nlist}")
    print(f"  Search probes (nprobe): {nprobe}")
    print(f"  PQ segments: 16, bits: 8")
    
    start_time = time.time()
    
    # Create index
    quantizer = faiss.IndexFlatIP(D)  # Inner product for cosine similarity
    index = faiss.IndexIVFPQ(quantizer, D, nlist, 16, 8)
    
    # Train
    print("  Training index...")
    train_start = time.time()
    index.train(embeddings.astype(np.float32))
    train_time = time.time() - train_start
    print(f"    Training time: {train_time:.1f}s")
    
    # Add embeddings
    print("  Adding embeddings...")
    add_start = time.time()
    index.add(embeddings.astype(np.float32))
    add_time = time.time() - add_start
    print(f"    Adding time: {add_time:.1f}s")
    
    # Set search parameters
    index.nprobe = nprobe
    
    build_time = time.time() - start_time
    
    print(f"\n✅ Index built!")
    print(f"  Total time: {build_time/60:.1f} minutes")
    print(f"  Index size: {index.ntotal:,} embeddings")
    
    return index


def save_index(index, positions, chr_names, output_dir: Path):
    """Save index and metadata to disk"""
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print(f"\nSaving index to disk...")
    
    # Save FAISS index
    index_file = output_dir / "genocache_v4_production.index"
    faiss.write_index(index, str(index_file))
    print(f"  ✅ FAISS index: {index_file}")
    
    # Save metadata
    metadata = {
        'positions': positions,
        'chr_names': chr_names,
        'config': {
            'window_size': 512,
            'stride': 32,
            'embed_dim': 128,
            'num_embeddings': len(positions)
        }
    }
    metadata_file = output_dir / "genocache_v4_production.metadata.pkl"
    with open(metadata_file, 'wb') as f:
        pickle.dump(metadata, f)
    print(f"  ✅ Metadata: {metadata_file}")
    
    # Calculate sizes
    index_size = index_file.stat().st_size / 1e9
    metadata_size = metadata_file.stat().st_size / 1e9
    total_size = index_size + metadata_size
    
    print(f"\n📊 Index Statistics:")
    print(f"  FAISS index: {index_size:.2f} GB")
    print(f"  Metadata: {metadata_size:.2f} GB")
    print(f"  Total: {total_size:.2f} GB")


def main():
    print("=" * 80)
    print("GenoCache V4 - Production Index Builder")
    print("=" * 80)
    print()
    
    # Configuration
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    CHECKPOINT_PATH = Path("/home/nebius/genocache/genocache-v4/models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt")
    OUTPUT_DIR = Path("/home/nebius/genocache/genocache-v4/indexes")
    
    # Primary chromosomes (chr1-22, X, Y)
    PRIMARY_CHROMOSOMES = [
        'NC_000001.11',  # chr1
        'NC_000002.12',  # chr2
        'NC_000003.12',  # chr3
        'NC_000004.12',  # chr4
        'NC_000005.10',  # chr5
        'NC_000006.12',  # chr6
        'NC_000007.14',  # chr7
        'NC_000008.11',  # chr8
        'NC_000009.12',  # chr9
        'NC_000010.11',  # chr10
        'NC_000011.10',  # chr11
        'NC_000012.12',  # chr12
        'NC_000013.11',  # chr13
        'NC_000014.9',   # chr14
        'NC_000015.10',  # chr15
        'NC_000016.10',  # chr16
        'NC_000017.11',  # chr17
        'NC_000018.10',  # chr18
        'NC_000019.10',  # chr19
        'NC_000020.11',  # chr20
        'NC_000021.9',   # chr21
        'NC_000022.11',  # chr22
        'NC_000023.11',  # chrX
        'NC_000024.10',  # chrY
    ]
    
    print("Configuration:")
    print(f"  Genome: {GENOME_PATH}")
    print(f"  Model: {CHECKPOINT_PATH.name}")
    print(f"  Chromosomes: {len(PRIMARY_CHROMOSOMES)} (chr1-22, X, Y)")
    print(f"  Output: {OUTPUT_DIR}")
    print()
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
        print(f"  Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print()
    
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
    
    checkpoint = torch.load(CHECKPOINT_PATH, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    num_params = sum(p.numel() for p in model.parameters())
    print(f"✅ Model loaded: {num_params:,} parameters")
    print()
    
    # Load genome
    genome = load_genome(GENOME_PATH, PRIMARY_CHROMOSOMES)
    
    if len(genome) < len(PRIMARY_CHROMOSOMES):
        print(f"⚠️  Warning: Only found {len(genome)} of {len(PRIMARY_CHROMOSOMES)} chromosomes")
        print(f"   Available: {list(genome.keys())}")
    
    total_bp = sum(len(seq) for seq in genome.values())
    print(f"✅ Loaded {len(genome)} chromosomes, {total_bp:,} bp total")
    print()
    
    # Encode genome
    embeddings, positions, chr_names = encode_genome(
        model, genome,
        window_size=512,
        stride=32,
        device=device,
        batch_size=512
    )
    
    # Build FAISS index
    index = build_faiss_index(embeddings, nprobe=64)
    
    # Save to disk
    save_index(index, positions, chr_names, OUTPUT_DIR)
    
    print("\n" + "=" * 80)
    print("✅ Production index build complete!")
    print("=" * 80)
    print()
    print("Index is ready for use in:")
    print("  1. Real data validation")
    print("  2. Head-to-head benchmarks vs minimap2")
    print("  3. Production deployment")


if __name__ == "__main__":
    main()
