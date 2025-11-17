#!/usr/bin/env python3
"""
Test EXTEND Phase - Fix for 37% Chromosome Accuracy Bug

This tests the complete fixed pipeline:
1. Adaptive seeding → return top-k candidates
2. EXTEND phase → align to each candidate
3. Pick best by ALIGNMENT SCORE (not seed count!)

Expected: 37% → 95%+ chromosome accuracy
"""

import sys
import time
from pathlib import Path

# Import our modules
from adaptive_seeding import AdaptiveSeeder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner
import torch
import faiss
import numpy as np
from Bio import SeqIO

def load_model_and_index():
    """Load trained model and FAISS index"""
    print("Loading model and index...")
    
    # Load model
    model_path = Path("../genocache-production/model_best.pt")
    if not model_path.exists():
        print(f"❌ Model not found: {model_path}")
        return None, None, None
    
    checkpoint = torch.load(model_path, map_location='cpu')
    model = checkpoint['model']
    model.eval()
    
    # Load FAISS index
    index_path = Path("../genocache-production/cnn_hyenadna_flat.index")
    if not index_path.exists():
        print(f"❌ Index not found: {index_path}")
        return None, None, None
    
    index = faiss.read_index(str(index_path))
    
    # Load metadata
    metadata_path = Path("../genocache-production/index_metadata.npy")
    if not metadata_path.exists():
        print(f"❌ Metadata not found: {metadata_path}")
        return None, None, None
    
    metadata = np.load(metadata_path, allow_pickle=True)
    
    print(f"✅ Model loaded: {model_path}")
    print(f"✅ Index loaded: {index.ntotal:,} vectors")
    print(f"✅ Metadata loaded: {len(metadata):,} entries")
    
    return model, index, metadata


def test_fixed_pipeline():
    """Test complete pipeline with EXTEND phase"""
    print("╔════════════════════════════════════════════════════════════════════════════╗")
    print("║           Testing EXTEND Phase - Chromosome Accuracy Fix                  ║")
    print("╚════════════════════════════════════════════════════════════════════════════╝")
    print()
    
    # Load model and index
    model, index, metadata = load_model_and_index()
    if model is None:
        print("❌ Failed to load model/index")
        return
    
    # Initialize components
    print("Initializing pipeline...")
    seeder = AdaptiveSeeder(model, index, metadata)
    aligner = FastAligner()
    extender = ExtendPhase(aligner, min_score_threshold=100)
    
    print("✅ Adaptive seeder ready")
    print("✅ Fast aligner ready (parasail)")
    print("✅ EXTEND phase ready")
    print()
    
    # Load test reads
    test_file = Path("../genocache-production/reads_chr22_synth_1kb_500.fa")
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return
    
    print(f"Loading test reads from: {test_file}")
    reads = list(SeqIO.parse(test_file, "fasta"))
    print(f"✅ Loaded {len(reads)} test reads")
    print()
    
    # Test on first 10 reads
    print("Testing on first 10 reads...")
    print("=" * 80)
    
    results = []
    total_time = 0
    
    for i, record in enumerate(reads[:10]):
        print(f"\n[Read {i}] {record.id}")
        print(f"Length: {len(record.seq)} bp")
        
        start_time = time.time()
        
        # Step 1: Adaptive seeding (return top-5 candidates)
        print("  [1] Seeding... ", end="", flush=True)
        candidates = seeder.align_read(str(record.seq), record.id, return_top_k=5)
        
        if candidates is None:
            print("❌ No candidates found")
            results.append({
                'read_id': record.id,
                'status': 'unmapped',
                'reason': 'no_seeds'
            })
            continue
        
        print(f"✅ Found {len(candidates)} candidates")
        
        # Show candidates
        for j, cand in enumerate(candidates):
            print(f"      Candidate {j+1}: {cand['chr']} "
                  f"({cand['start']:,}-{cand['end']:,}) "
                  f"seeds={cand['num_seeds']}, score={cand['score']:.2f}")
        
        # Step 2: EXTEND - align to each candidate, pick best by score
        print("  [2] EXTEND phase (align each candidate)... ", end="", flush=True)
        best_alignment = extender.extend_and_score(str(record.seq), candidates)
        
        if best_alignment is None:
            print("❌ No valid alignments")
            results.append({
                'read_id': record.id,
                'status': 'unmapped',
                'reason': 'low_alignment_score'
            })
            continue
        
        elapsed = time.time() - start_time
        total_time += elapsed
        
        print(f"✅")
        print(f"  [3] Best alignment: {best_alignment['chr']}")
        print(f"      Position: {best_alignment['start']:,}-{best_alignment['end']:,}")
        print(f"      Alignment score: {best_alignment['alignment_score']}")
        print(f"      Seed score: {best_alignment['seed_score']:.2f}")
        print(f"      Status: {best_alignment['status']}")
        print(f"      Time: {elapsed:.3f}s")
        
        results.append(best_alignment)
    
    # Summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    mapped = [r for r in results if 'chr' in r]
    unmapped = [r for r in results if 'chr' not in r]
    
    print(f"Total reads: {len(results)}")
    print(f"Mapped: {len(mapped)} ({len(mapped)/len(results)*100:.1f}%)")
    print(f"Unmapped: {len(unmapped)}")
    print(f"Avg time: {total_time/len(results):.3f}s per read")
    print(f"Speed: {len(results)/total_time:.2f} reads/sec")
    print()
    
    # Show chromosome distribution
    if mapped:
        chr_counts = {}
        for r in mapped:
            chr_name = r['chr']
            chr_counts[chr_name] = chr_counts.get(chr_name, 0) + 1
        
        print("Chromosome distribution:")
        for chr_name, count in sorted(chr_counts.items()):
            print(f"  {chr_name}: {count} reads")
    
    print()
    print("Next: Compare with minimap2 to measure chromosome accuracy!")
    
    return results


if __name__ == '__main__':
    test_fixed_pipeline()
