#!/usr/bin/env python3
"""
Phase 5.3: Build FAISS index from full genome vectors
Creates optimized index for fast similarity search
"""

import numpy as np
import faiss
import os
import time
from glob import glob

# Configuration
CHECKPOINT_DIR = "/home/nebius/genocache/encoding_checkpoints_full"
OUTPUT_INDEX = "/home/nebius/genocache/grch38_full_index.faiss"
OUTPUT_POSITIONS = "/home/nebius/genocache/grch38_full_positions.npy"
D = 256  # Embedding dimension
USE_GPU = True

print("="*70)
print("Phase 5.3: Build FAISS Index from Full Genome")
print("="*70)

# Step 1: Load all checkpoint files
print("\n[1/5] Loading encoded vectors from checkpoints...")

checkpoint_files = sorted(glob(f"{CHECKPOINT_DIR}/chr_*.npz"))
print(f"  Found {len(checkpoint_files)} checkpoint files")

if len(checkpoint_files) == 0:
    print("ERROR: No checkpoint files found!")
    print(f"  Looking in: {CHECKPOINT_DIR}")
    exit(1)

all_vectors = []
all_positions = []

for i, ckpt_file in enumerate(checkpoint_files, 1):
    print(f"  Loading {i}/{len(checkpoint_files)}: {os.path.basename(ckpt_file)}")
    data = np.load(ckpt_file)
    
    vectors = data['vectors']
    positions = data['positions']
    
    print(f"    Vectors: {vectors.shape}, Positions: {positions.shape}")
    
    all_vectors.append(vectors)
    all_positions.append(positions)

print("\n[2/5] Concatenating all vectors...")
vectors = np.vstack(all_vectors)
positions = np.concatenate(all_positions)

print(f"  ✓ Total vectors: {len(vectors):,}")
print(f"  ✓ Vector shape: {vectors.shape}")
print(f"  ✓ Memory: {vectors.nbytes / 1e9:.2f} GB")

# Save positions
print(f"\n[3/5] Saving position map...")
np.save(OUTPUT_POSITIONS, positions)
print(f"  ✓ Saved: {OUTPUT_POSITIONS}")
print(f"  ✓ Size: {os.path.getsize(OUTPUT_POSITIONS) / 1e6:.1f} MB")

# Step 2: Build FAISS index
print("\n[4/5] Building FAISS index...")

n_vectors = len(vectors)

# Choose index type based on dataset size
if n_vectors < 1_000_000:
    # Small: Use Flat index (exact search)
    print("  Using Flat index (exact search)")
    index = faiss.IndexFlatL2(D)
    
elif n_vectors < 10_000_000:
    # Medium: Use IVF with clustering
    print("  Using IVF index with clustering")
    nlist = min(int(np.sqrt(n_vectors)), 50000)  # Number of clusters
    quantizer = faiss.IndexFlatL2(D)
    index = faiss.IndexIVFFlat(quantizer, D, nlist)
    
    print(f"  Training with {nlist:,} clusters...")
    # Sample for training if too large
    if n_vectors > 1_000_000:
        train_idx = np.random.choice(n_vectors, 1_000_000, replace=False)
        train_vectors = vectors[train_idx].astype(np.float32)
    else:
        train_vectors = vectors.astype(np.float32)
    
    index.train(train_vectors)
    print("  ✓ Training complete")
    
else:
    # Large: Use IVF with Product Quantization (compression)
    print("  Using IVF-PQ index with compression")
    nlist = min(int(np.sqrt(n_vectors)), 100000)
    m = 64  # Number of subquantizers
    nbits = 8  # Bits per subquantizer
    
    quantizer = faiss.IndexFlatL2(D)
    index = faiss.IndexIVFPQ(quantizer, D, nlist, m, nbits)
    
    print(f"  Training with {nlist:,} clusters, PQ: m={m}, nbits={nbits}...")
    # Sample for training
    train_idx = np.random.choice(n_vectors, min(2_000_000, n_vectors), replace=False)
    train_vectors = vectors[train_idx].astype(np.float32)
    
    index.train(train_vectors)
    print("  ✓ Training complete")

# Add vectors to index
print(f"\n  Adding {n_vectors:,} vectors to index...")
start_time = time.time()

# Add in batches to avoid memory issues
batch_size = 100_000
for i in range(0, n_vectors, batch_size):
    end_idx = min(i + batch_size, n_vectors)
    batch = vectors[i:end_idx].astype(np.float32)
    index.add(batch)
    
    if (i // batch_size) % 10 == 0:
        progress = (end_idx / n_vectors) * 100
        elapsed = time.time() - start_time
        rate = end_idx / elapsed if elapsed > 0 else 0
        print(f"    Progress: {progress:.1f}% ({end_idx:,}/{n_vectors:,}) @ {rate:,.0f} vectors/sec")

elapsed = time.time() - start_time
print(f"  ✓ Added all vectors in {elapsed:.1f}s ({n_vectors/elapsed:,.0f} vectors/sec)")

# Set search parameters
if hasattr(index, 'nprobe'):
    index.nprobe = 64  # Number of clusters to search
    print(f"  ✓ Set nprobe = {index.nprobe}")

# Step 3: Save index
print(f"\n[5/5] Saving FAISS index...")
faiss.write_index(index, OUTPUT_INDEX)
print(f"  ✓ Saved: {OUTPUT_INDEX}")

index_size = os.path.getsize(OUTPUT_INDEX)
print(f"  ✓ Size: {index_size / 1e6:.1f} MB")
print(f"  ✓ Compression ratio: {vectors.nbytes / index_size:.1f}x")

# Test search
print("\n[Test] Running test search...")
test_query = vectors[:5].astype(np.float32)
k = 10
distances, indices = index.search(test_query, k)

print(f"  ✓ Test search successful")
print(f"  Query shape: {test_query.shape}")
print(f"  Result distances: {distances[0, :3]}")
print(f"  Result indices: {indices[0, :3]}")

print("\n" + "="*70)
print("✓ Phase 5.3 Complete!")
print("="*70)
print(f"\nIndex Statistics:")
print(f"  Total vectors:     {n_vectors:,}")
print(f"  Dimension:         {D}")
print(f"  Index size:        {index_size / 1e6:.1f} MB")
print(f"  Vectors size:      {vectors.nbytes / 1e9:.2f} GB")
print(f"  Compression:       {vectors.nbytes / index_size:.1f}x")
print(f"  Index type:        {type(index).__name__}")

print(f"\nOutput files:")
print(f"  Index:     {OUTPUT_INDEX}")
print(f"  Positions: {OUTPUT_POSITIONS}")

print("\nNext steps:")
print("  1. Test alignment with full genome index")
print("  2. Phase 5.4: Validate on HG002 dataset")
print("  3. Phase 5.5: Compare with minimap2")
