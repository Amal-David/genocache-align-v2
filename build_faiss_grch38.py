#!/usr/bin/env python3
"""
Build full-genome FAISS index for GRCh38
"""
import os, sys, math, time
import numpy as np
import faiss

REF_VECTORS = "grch38_vectors.npy"
REF_POSITIONS = "grch38_positions.npy"
REF_CHROMOSOMES = "grch38_chromosomes.npy"
OUT_INDEX = "faiss_grch38.idx"

# FAISS parameters for large index
MIN_NLIST = 128
NPROBE = 64  # higher for better recall on large index
PQ_M = 32  # 256 / 8 = 32 subquantizers

def main():
    print("=" * 70)
    print("Building Full-Genome FAISS Index (GRCh38)")
    print("=" * 70)
    
    if not os.path.exists(REF_VECTORS):
        print(f"\nERROR: {REF_VECTORS} not found")
        print("Please run encode_grch38.py first")
        sys.exit(1)
    
    # Load vectors
    print(f"\n[1/5] Loading reference vectors...")
    start_time = time.time()
    ref_vecs = np.load(REF_VECTORS, mmap_mode='r')  # memory-mapped for large files
    N, DIM = ref_vecs.shape
    print(f"  ✓ Shape: {N:,} vectors × {DIM} dimensions")
    print(f"  ✓ Size: {ref_vecs.nbytes / 1024**3:.2f} GB")
    
    # Load metadata
    print(f"\n[2/5] Loading position metadata...")
    ref_pos = np.load(REF_POSITIONS)
    ref_chr = np.load(REF_CHROMOSOMES)
    print(f"  ✓ Positions: {len(ref_pos):,}")
    print(f"  ✓ Chromosomes: {len(np.unique(ref_chr))} unique")
    
    # Determine parameters
    nlist = max(MIN_NLIST, int(math.sqrt(N) * 2))
    nlist = min(nlist, N // 100)
    
    print(f"\n[3/5] Index configuration:")
    print(f"  Index type: IVF-PQ")
    print(f"  nlist (clusters): {nlist:,}")
    print(f"  nprobe (search): {NPROBE}")
    print(f"  PQ subquantizers: {PQ_M}")
    print(f"  Metric: Inner Product (cosine similarity)")
    
    # Build index
    print(f"\n[4/5] Building index...")
    
    # Create quantizer
    print(f"  Creating quantizer...")
    quantizer = faiss.IndexFlatIP(DIM)
    index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)
    index.metric_type = faiss.METRIC_INNER_PRODUCT
    
    # Train
    train_size = min(500_000, N)  # larger training set
    train_size = max(train_size, 20 * nlist)
    if train_size > N:
        train_size = N
    
    print(f"  Training with {train_size:,} samples...")
    rng = np.random.RandomState(42)
    train_indices = rng.permutation(N)[:train_size]
    
    # Load training data in chunks if needed
    if N > 10_000_000:  # if very large, load in chunks
        train_data = ref_vecs[train_indices].copy()
    else:
        train_data = ref_vecs[train_indices]
    
    train_start = time.time()
    index.train(train_data.astype('float32'))
    train_time = time.time() - train_start
    print(f"  ✓ Training completed in {train_time:.1f}s")
    
    # Add vectors in batches
    print(f"  Adding {N:,} vectors...")
    add_start = time.time()
    
    BATCH_SIZE = 100_000
    for i in range(0, N, BATCH_SIZE):
        end_i = min(i + BATCH_SIZE, N)
        batch = ref_vecs[i:end_i]
        if isinstance(batch, np.memmap):
            batch = np.array(batch, dtype='float32')
        index.add(batch.astype('float32'))
        
        if (i + BATCH_SIZE) % 1_000_000 == 0 or end_i == N:
            progress = end_i / N * 100
            elapsed = time.time() - add_start
            rate = end_i / elapsed
            eta = (N - end_i) / rate if rate > 0 else 0
            print(f"    Progress: {end_i:,}/{N:,} ({progress:.1f}%) - {rate:.0f} vec/s - ETA: {eta:.0f}s")
    
    add_time = time.time() - add_start
    print(f"  ✓ Added all vectors in {add_time:.1f}s ({N/add_time:.0f} vec/s)")
    
    # Set search parameters
    index.nprobe = NPROBE
    
    # Save
    print(f"\n[5/5] Saving index...")
    save_start = time.time()
    faiss.write_index(index, OUT_INDEX)
    save_time = time.time() - save_start
    
    file_size = os.path.getsize(OUT_INDEX) / (1024**3)
    print(f"  ✓ Saved: {OUT_INDEX}")
    print(f"  ✓ Size: {file_size:.2f} GB")
    print(f"  ✓ Save time: {save_time:.1f}s")
    
    # Test search
    print(f"\n[Test] Running sample searches...")
    test_queries = 100
    test_indices = rng.choice(N, size=test_queries, replace=False)
    test_vecs = ref_vecs[test_indices]
    if isinstance(test_vecs, np.memmap):
        test_vecs = np.array(test_vecs, dtype='float32')
    
    search_start = time.time()
    D, I = index.search(test_vecs.astype('float32'), 10)
    search_time = time.time() - search_start
    
    # Check recall (top-1 should match self)
    recall_1 = np.sum(I[:, 0] == test_indices) / test_queries
    
    print(f"  ✓ Search time: {search_time*1000:.1f}ms for {test_queries} queries")
    print(f"  ✓ Throughput: {test_queries/search_time:.0f} queries/sec")
    print(f"  ✓ Self-recall@1: {recall_1*100:.1f}% (should be ~100%)")
    
    total_time = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"Index building complete!")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Index contains {N:,} vectors from {len(np.unique(ref_chr))} chromosomes")
    print("=" * 70)

if __name__ == "__main__":
    main()
