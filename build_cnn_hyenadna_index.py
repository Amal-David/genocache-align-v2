#!/usr/bin/env python3
"""
Build FAISS Flat Index for CNN-HyenaDNA vectors
Exact search - no quantization loss
"""

import numpy as np
import faiss
import time
from pathlib import Path

# Config
VECTORS_FILE = "cnn_hyenadna_vectors.npy"
INDEX_OUTPUT = "cnn_hyenadna_flat.index"
DIMENSION = 256

print("="*80)
print("CNN-HYENADNA FAISS FLAT INDEX BUILDER")
print("="*80)
print()

# Load vectors
print(f"Loading vectors from {VECTORS_FILE}...")
start_time = time.time()
vectors = np.load(VECTORS_FILE)
load_time = time.time() - start_time

print(f"✓ Loaded vectors in {load_time:.1f} seconds")
print(f"  Shape: {vectors.shape}")
print(f"  Total vectors: {len(vectors):,}")
print(f"  Dimension: {vectors.shape[1]}")
print(f"  Size: {vectors.nbytes / 1e9:.2f} GB")
print(f"  Dtype: {vectors.dtype}")
print()

# Verify dimension
assert vectors.shape[1] == DIMENSION, f"Expected dimension {DIMENSION}, got {vectors.shape[1]}"

# Convert to float32 if needed
if vectors.dtype != np.float32:
    print(f"Converting from {vectors.dtype} to float32...")
    vectors = vectors.astype(np.float32)
    print("✓ Converted")
    print()

# Normalize vectors for cosine similarity
print("Normalizing vectors (L2 norm)...")
start_time = time.time()
faiss.normalize_L2(vectors)
normalize_time = time.time() - start_time
print(f"✓ Normalized in {normalize_time:.1f} seconds")
print()

# Build Flat index (exact search)
print(f"Building Flat index (exact search, no quantization)...")
print(f"  Index type: IndexFlatIP (Inner Product)")
print(f"  Vectors: {len(vectors):,}")
print(f"  Dimension: {DIMENSION}")
print()

start_time = time.time()

# Create flat index for inner product (cosine similarity after normalization)
index = faiss.IndexFlatIP(DIMENSION)

# Add vectors
print("Adding vectors to index...")
index.add(vectors)
add_time = time.time() - start_time

print(f"✓ Index built in {add_time:.1f} seconds")
print(f"  Total vectors in index: {index.ntotal:,}")
print()

# Verify index
print("Verifying index...")
assert index.ntotal == len(vectors), "Vector count mismatch!"
print("✓ Index verified")
print()

# Save index
print(f"Saving index to {INDEX_OUTPUT}...")
start_time = time.time()
faiss.write_index(index, INDEX_OUTPUT)
save_time = time.time() - start_time
print(f"✓ Saved in {save_time:.1f} seconds")
print()

# Get file size
index_size = Path(INDEX_OUTPUT).stat().st_size / 1e9
print(f"Index file size: {index_size:.2f} GB")
print()

# Test search
print("Testing index with sample query...")
test_vector = vectors[0:1].copy()
faiss.normalize_L2(test_vector)
distances, indices = index.search(test_vector, k=5)
print(f"✓ Search works! Top 5 results:")
for i, (dist, idx) in enumerate(zip(distances[0], indices[0])):
    print(f"  {i+1}. Index {idx:,} (distance: {dist:.4f})")
print()

print("="*80)
print("✅ INDEX BUILD COMPLETE!")
print("="*80)
print()
print(f"Summary:")
print(f"  Vectors: {len(vectors):,}")
print(f"  Dimension: {DIMENSION}")
print(f"  Index type: Flat (exact search)")
print(f"  Index size: {index_size:.2f} GB")
print(f"  Output file: {INDEX_OUTPUT}")
print()
print(f"Total time: {load_time + normalize_time + add_time + save_time:.1f} seconds")
print()
print("Next steps:")
print("  1. Test search accuracy with validation queries")
print("  2. Integrate into 6-model ensemble")
print("  3. Run full genome validation")
print()
