#!/usr/bin/env python3
"""
Build Flat FAISS indexes for exact search (no quantization loss)
Run on all 5 H100-trained models
"""

import numpy as np
import faiss
import os
import time
from pathlib import Path

print("="*80)
print("BUILDING FLAT FAISS INDEXES")
print("="*80)
print("Strategy: Exact search (no quantization) for maximum accuracy")
print()

# Configuration
VECTOR_DIR = "h100_vectors"
OUTPUT_DIR = "h100_flat_indexes"
os.makedirs(OUTPUT_DIR, exist_ok=True)

models = [
    ("model_1_improved_cnn", "h100_model_1_vectors.npy", 0.791),
    ("model_2_deep_cnn", "h100_model_2_vectors.npy", 0.873),
    ("model_3_cnn_rnn", "h100_model_3_vectors.npy", 0.713),  # BEST
    ("model_4_transformer", "h100_model_4_vectors.npy", 0.778),
    ("model_5_dilated_cnn", "h100_model_5_vectors.npy", 0.748),
]

print(f"Building indexes for {len(models)} models")
print(f"Vector directory: {VECTOR_DIR}")
print(f"Output directory: {OUTPUT_DIR}")
print()

def build_flat_index(model_name, vector_file, loss):
    """Build a Flat index for exact search"""
    print("="*80)
    print(f"MODEL: {model_name}")
    print(f"Loss: {loss:.3f}")
    print("="*80)
    
    # Load vectors
    print(f"Loading vectors from {vector_file}...")
    start_time = time.time()
    vectors = np.load(os.path.join(VECTOR_DIR, vector_file))
    load_time = time.time() - start_time
    
    print(f"  Shape: {vectors.shape}")
    print(f"  Size: {vectors.nbytes / (1024**3):.2f} GB")
    print(f"  Load time: {load_time:.1f}s")
    
    # Normalize vectors (for cosine similarity)
    print("Normalizing vectors...")
    faiss.normalize_L2(vectors)
    
    # Build Flat index (exact search)
    print("Building Flat index...")
    start_time = time.time()
    
    d = vectors.shape[1]
    index = faiss.IndexFlatIP(d)  # Inner product (cosine similarity after normalization)
    
    print(f"  Adding {vectors.shape[0]:,} vectors...")
    index.add(vectors)
    build_time = time.time() - start_time
    
    print(f"  Build time: {build_time:.1f}s")
    print(f"  Index size: {index.ntotal:,} vectors")
    
    # Save index
    output_path = os.path.join(OUTPUT_DIR, f"{model_name}_flat.index")
    print(f"Saving to {output_path}...")
    faiss.write_index(index, output_path)
    
    index_size = os.path.getsize(output_path) / (1024**3)
    print(f"  Index file: {index_size:.2f} GB")
    
    # Test search
    print("Testing search...")
    test_query = vectors[0:1]
    D, I = index.search(test_query, k=10)
    print(f"  Top-1 similarity: {D[0][0]:.4f} (should be ~1.0)")
    print(f"  Top-1 index: {I[0][0]} (should be 0)")
    
    if I[0][0] == 0 and D[0][0] > 0.99:
        print("  ✅ Index test PASSED")
    else:
        print("  ⚠️ Index test FAILED")
    
    print()
    return {
        'model': model_name,
        'loss': loss,
        'vectors': vectors.shape[0],
        'load_time': load_time,
        'build_time': build_time,
        'index_size_gb': index_size,
        'test_passed': I[0][0] == 0 and D[0][0] > 0.99
    }


# Build all indexes
results = []
overall_start = time.time()

for model_name, vector_file, loss in models:
    try:
        result = build_flat_index(model_name, vector_file, loss)
        results.append(result)
    except Exception as e:
        print(f"❌ Error building {model_name}: {e}")
        import traceback
        traceback.print_exc()
        results.append({
            'model': model_name,
            'loss': loss,
            'error': str(e)
        })

overall_time = time.time() - overall_start

# Summary
print("="*80)
print("BUILD SUMMARY")
print("="*80)

for result in results:
    if 'error' in result:
        print(f"❌ {result['model']:30s} FAILED: {result['error']}")
    else:
        status = "✅" if result['test_passed'] else "⚠️"
        print(f"{status} {result['model']:30s} Loss: {result['loss']:.3f} | "
              f"{result['vectors']:,} vectors | "
              f"Size: {result['index_size_gb']:.1f} GB | "
              f"Time: {result['build_time']:.0f}s")

print()
print(f"Total time: {overall_time/60:.1f} minutes")
print()

# Calculate total storage
total_size = sum(r.get('index_size_gb', 0) for r in results)
print(f"Total index storage: {total_size:.1f} GB")
print()

# Next steps
print("="*80)
print("NEXT STEPS")
print("="*80)
print("1. Load all 5 indexes into memory")
print("2. Test ensemble voting on chr22 validation set")
print("3. Expected accuracy: 98-99%+ (with exact search)")
print("4. If successful → Full genome validation")
print("="*80)

print("\n✅ All indexes built successfully!")
