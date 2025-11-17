#!/usr/bin/env python3
"""
Resume encoding from checkpoints and assemble final files
"""
import os, sys, glob
import numpy as np

CHECKPOINT_DIR = "encoding_checkpoints"
OUT_VECTORS = "grch38_vectors.npy"
OUT_POSITIONS = "grch38_positions.npy"
OUT_CHROMOSOMES = "grch38_chromosomes.npy"

def main():
    print("=" * 70)
    print("Assembling GRCh38 from Checkpoints")
    print("=" * 70)
    
    checkpoint_files = sorted(glob.glob(f"{CHECKPOINT_DIR}/chr_*.npz"))
    
    if not checkpoint_files:
        print(f"\nERROR: No checkpoint files found in {CHECKPOINT_DIR}/")
        print("Run encode_grch38_robust.py first")
        sys.exit(1)
    
    print(f"\nFound {len(checkpoint_files)} checkpoint files")
    
    # Load all checkpoints
    print(f"\nLoading checkpoints...")
    all_vectors = []
    all_positions = []
    all_chromosomes = []
    
    for i, ckpt in enumerate(checkpoint_files):
        chr_name = os.path.basename(ckpt)
        print(f"  [{i+1}/{len(checkpoint_files)}] {chr_name}")
        
        try:
            data = np.load(ckpt)
            all_vectors.append(data['vectors'])
            all_positions.extend(data['positions'].tolist())
            all_chromosomes.extend(data['chromosomes'].tolist())
            
            print(f"      ✓ {data['vectors'].shape[0]:,} vectors")
        except Exception as e:
            print(f"      ✗ Failed to load: {e}")
            continue
    
    print(f"\n  Loaded {len(all_vectors)} chromosomes")
    
    # Concatenate
    print(f"\nConcatenating arrays...")
    try:
        vectors = np.vstack(all_vectors).astype('float32')
        positions = np.array(all_positions, dtype=np.int64)
        chromosomes = np.array(all_chromosomes, dtype=np.int32)
        
        print(f"  ✓ Total vectors: {len(vectors):,}")
        print(f"  ✓ Vector shape: {vectors.shape}")
        print(f"  ✓ Unique chromosomes: {len(np.unique(chromosomes))}")
    except Exception as e:
        print(f"  ✗ Concatenation failed: {e}")
        sys.exit(1)
    
    # Save final files
    print(f"\nSaving final files...")
    
    try:
        print(f"  Saving {OUT_VECTORS}...")
        np.save(OUT_VECTORS, vectors)
        size_gb = os.path.getsize(OUT_VECTORS) / 1024**3
        print(f"    ✓ Saved ({size_gb:.2f} GB)")
        
        print(f"  Saving {OUT_POSITIONS}...")
        np.save(OUT_POSITIONS, positions)
        size_mb = os.path.getsize(OUT_POSITIONS) / 1024**2
        print(f"    ✓ Saved ({size_mb:.1f} MB)")
        
        print(f"  Saving {OUT_CHROMOSOMES}...")
        np.save(OUT_CHROMOSOMES, chromosomes)
        size_mb = os.path.getsize(OUT_CHROMOSOMES) / 1024**2
        print(f"    ✓ Saved ({size_mb:.1f} MB)")
        
        print(f"\n{'='*70}")
        print("Assembly complete!")
        print(f"Total vectors: {len(vectors):,}")
        print("=" * 70)
        
    except Exception as e:
        print(f"  ✗ Save failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
