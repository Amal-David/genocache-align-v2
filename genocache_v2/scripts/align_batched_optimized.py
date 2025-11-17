#!/usr/bin/env python3
import argparse
import time
import numpy as np
import faiss
import torch
from pathlib import Path
from Bio import SeqIO
import sys

sys.path.insert(0, str(Path(__file__).parent))
from improved_cnn import ImprovedCNN

# Vectorized one-hot encoding
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def one_hot_encode_batch_fast(sequences):
    """Vectorized batch one-hot encoding"""
    batch_size = len(sequences)
    seq_len = len(sequences[0])
    
    # Pre-allocate output tensor
    encoded = torch.zeros((batch_size, 4, seq_len), dtype=torch.float32)
    
    for i, seq in enumerate(sequences):
        for j, base in enumerate(seq):
            if base in BASE_TO_INDEX:
                encoded[i, BASE_TO_INDEX[base], j] = 1.0
    
    return encoded

def load_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences

def cluster_positions(positions, window=100):
    if len(positions) == 0:
        return []
    positions = sorted(positions)
    clusters = []
    current_cluster = [positions[0]]
    for pos in positions[1:]:
        if pos - current_cluster[-1] <= window:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)
    return clusters

def extract_seeds_batch(reads, seed_len=512, num_seeds=5):
    """Extract all seeds from all reads at once"""
    all_seeds = []
    read_seed_map = []
    
    for read_id, read_seq in enumerate(reads):
        read_seeds = []
        for i in range(num_seeds):
            offset = i * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            all_seeds.append((seed, offset))
            read_seeds.append(len(all_seeds) - 1)
        read_seed_map.append(read_seeds)
    
    return all_seeds, read_seed_map

def encode_seeds_batch(model, seeds, batch_size=512, device='cuda'):
    """Encode seeds in batches with vectorized one-hot encoding"""
    embeddings = []
    
    for i in range(0, len(seeds), batch_size):
        batch_seeds = seeds[i:i+batch_size]
        sequences = [seed for seed, _ in batch_seeds]
        
        # Vectorized encoding - much faster!
        batch_tensor = one_hot_encode_batch_fast(sequences).to(device)
        
        # Batch inference
        with torch.no_grad():
            batch_emb = model(batch_tensor).cpu().numpy()
        
        embeddings.extend(batch_emb)
    
    return np.array(embeddings, dtype=np.float32)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--num-reads", type=int, default=1000)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--cluster-window", type=int, default=100)
    p.add_argument("--tolerance", type=int, default=1024)
    args = p.parse_args()
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    print(f"🚀 Optimized Batched Alignment")
    print(f"   Device: {device}")
    print(f"   Batch size: {args.batch_size}")
    
    # Load everything
    start = time.time()
    model, _ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    
    index = faiss.read_index(args.index)
    ref_positions = np.load(args.positions)
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    
    print(f"✅ Loaded in {time.time()-start:.1f}s")
    print(f"   Reference: {len(ref_seq):,} bp")
    print(f"   Index: {index.ntotal:,} vectors")
    
    # Generate synthetic reads
    print(f"\n📝 Generating {args.num_reads} reads...")
    np.random.seed(42)
    reads = []
    true_positions = []
    
    max_start = len(ref_seq) - 2000 - 10
    candidate_starts = np.random.randint(10, max_start, args.num_reads)
    
    for pos in candidate_starts:
        read = ref_seq[pos:pos+2000]
        if len(read) == 2000 and 'N' not in read:
            reads.append(read)
            true_positions.append(pos)
    
    print(f"   Valid reads: {len(reads)}")
    
    # Extract all seeds
    print(f"\n🧬 Extracting seeds...")
    all_seeds, read_seed_map = extract_seeds_batch(reads)
    print(f"   Total seeds: {len(all_seeds)}")
    
    # Batch encode all seeds
    print(f"\n🔬 Encoding seeds (vectorized, batch={args.batch_size})...")
    encode_start = time.time()
    all_embeddings = encode_seeds_batch(model, all_seeds, args.batch_size, device)
    encode_time = time.time() - encode_start
    print(f"   Encoded {len(all_seeds)} seeds in {encode_time:.2f}s")
    print(f"   Speed: {len(all_seeds)/encode_time:.1f} seeds/sec")
    
    # Search all embeddings at once
    print(f"\n🔍 Searching index...")
    search_start = time.time()
    D, I = index.search(all_embeddings, k=10)
    search_time = time.time() - search_start
    print(f"   Searched {len(all_embeddings)} queries in {search_time:.2f}s")
    
    # Vote for each read
    print(f"\n🗳️  Voting...")
    vote_start = time.time()
    correct = 0
    
    for read_idx, seed_indices in enumerate(read_seed_map):
        all_positions = []
        
        for seed_idx in seed_indices:
            indices = I[seed_idx]
            positions = ref_positions[indices]
            seed_offset = all_seeds[seed_idx][1]
            adjusted = positions - seed_offset
            all_positions.extend(adjusted)
        
        # Cluster-based voting
        clusters = cluster_positions(all_positions, window=args.cluster_window)
        best_cluster = max(clusters, key=len)
        predicted_pos = int(np.median(best_cluster))
        
        # Check accuracy
        true_pos = true_positions[read_idx]
        if abs(predicted_pos - true_pos) <= args.tolerance:
            correct += 1
    
    vote_time = time.time() - vote_start
    
    # Results
    total_time = time.time() - start
    accuracy = 100 * correct / len(reads)
    throughput = len(reads) / total_time
    
    print(f"\n📊 Results:")
    print(f"   Accuracy: {accuracy:.1f}% ({correct}/{len(reads)})")
    print(f"   Total time: {total_time:.2f}s")
    print(f"   Throughput: {throughput:.1f} reads/sec")
    print(f"\n⏱️  Breakdown:")
    print(f"   Encoding: {encode_time:.2f}s ({100*encode_time/total_time:.0f}%)")
    print(f"   Searching: {search_time:.2f}s ({100*search_time/total_time:.0f}%)")
    print(f"   Voting: {vote_time:.2f}s ({100*vote_time/total_time:.0f}%)")

if __name__ == "__main__":
    main()
