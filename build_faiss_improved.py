#!/usr/bin/env python3
"""
Build optimized FAISS index with improved vectors
"""
import os, sys, math
import numpy as np
import faiss

REF_VECTORS = "ref_vectors_improved.npy"
REF_POSITIONS = "ref_positions_improved.npy"
OUT_INDEX = "faiss_index_improved.idx"

# FAISS parameters
MIN_NLIST = 64
NPROBE = 32  # increased for better recall
PQ_M = 32  # 256 / 8 = 32 subquantizers

def main():
    print("=" * 60)
    print("Building Optimized FAISS Index")
    print("=" * 60)
    
    if not os.path.exists(REF_VECTORS):
        print(f"ERROR: {REF_VECTORS} not found")
        sys.exit(1)
    
    print(f"\n[1/4] Loading reference vectors...")
    ref_vecs = np.load(REF_VECTORS).astype("float32")
    N, DIM = ref_vecs.shape
    print(f"  ✓ Loaded {N:,} vectors of dimension {DIM}")
    
    # Determine index parameters
    nlist = max(MIN_NLIST, int(math.sqrt(N) * 2))
    nlist = min(nlist, N // 100)  # at least 100 vectors per cluster
    
    print(f"\n[2/4] Index parameters:")
    print(f"  nlist (clusters): {nlist}")
    print(f"  nprobe (search): {NPROBE}")
    print(f"  PQ_M (subquantizers): {PQ_M}")
    
    # Build index
    print(f"\n[3/4] Building IVF-PQ index...")
    quantizer = faiss.IndexFlatIP(DIM)
    index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)
    index.metric_type = faiss.METRIC_INNER_PRODUCT
    
    # Train
    train_size = min(200_000, N)
    train_size = max(train_size, 10 * nlist)
    if train_size > N:
        train_size = N
    
    print(f"  Training with {train_size:,} samples...")
    rnd = np.random.RandomState(42)
    perm = rnd.permutation(N)[:train_size]
    index.train(ref_vecs[perm])
    
    # Add vectors
    print(f"  Adding {N:,} vectors...")
    index.add(ref_vecs)
    
    # Set nprobe
    index.nprobe = NPROBE
    
    print(f"\n[4/4] Saving index...")
    faiss.write_index(index, OUT_INDEX)
    
    file_size = os.path.getsize(OUT_INDEX) / (1024**2)
    print(f"  ✓ Saved: {OUT_INDEX} ({file_size:.1f} MB)")
    
    # Test search
    print(f"\n[Test] Running sample search...")
    test_vec = ref_vecs[0:1]
    D, I = index.search(test_vec, 10)
    print(f"  ✓ Top-1 match: index={I[0,0]} score={D[0,0]:.4f}")
    
    print("\n" + "=" * 60)
    print("Index building complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
