#!/usr/bin/env python3
"""
Test Chr22 NAL model vs minimap2 on synthetic reads with ground truth
"""

import sys
import torch
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from seeding_chr22 import Chr22Seeding
from chaining_nal import NALChaining

def load_reads(fastq_path):
    """Load FASTQ reads"""
    reads = []
    with open(fastq_path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # quality
            
            read_id = header.strip()[1:].split()[0]
            reads.append({'id': read_id, 'seq': seq})
    
    return reads

def load_ground_truth(csv_path):
    """Load ground truth positions"""
    truth = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            truth[row['read_id']] = {
                'chr': row['chr'],
                'start': int(row['start']),
                'end': int(row['end'])
            }
    return truth

def main():
    print("="*80)
    print("CHR22 NAL TESTING - SYNTHETIC READS WITH GROUND TRUTH")
    print("="*80)
    print()
    
    # Initialize NAL seeder and chainer
    print("Loading NAL components...")
    seeder = Chr22Seeding(
        model_path='models/chr22_nal_512bp_final.pt',
        index_path='indexes/chr22_nal_512bp_stride32.index',
        positions_path='indexes/chr22_nal_512bp_stride32.positions.npz',
        seed_len=512,
        K=32,
        device='cuda'
    )
    
    chainer = NALChaining(tolerance=1000, min_chain_score=3, top_k=5)
    
    # Load test data
    print("\nLoading test data...")
    reads = load_reads('data/chr22_test_reads.fastq')
    ground_truth = load_ground_truth('data/chr22_ground_truth.csv')
    print(f"  Loaded {len(reads)} synthetic reads")
    print(f"  Ground truth: {len(ground_truth)} positions")
    
    # Align reads
    print("\nAligning reads with NAL...")
    results = []
    
    for i, read in enumerate(reads):
        if (i + 1) % 50 == 0:
            print(f"  Progress: {i+1}/{len(reads)}")
        
        # Get anchors (6 seeds, K=32 neighbors per seed)
        anchors = seeder.get_anchors(read['seq'], num_seeds=6, K=32)
        
        if not anchors:
            results.append({
                'read_id': read['id'],
                'mapped': False,
                'chr': None,
                'pos': None
            })
            continue
        
        # Chain anchors
        chains, _ = chainer.chain_with_rescue_check(
            anchors,
            read_len=len(read['seq']),
            num_seeds=6,
            seed_len=512
        )
        
        if chains:
            best_chain = chains[0]
            results.append({
                'read_id': read['id'],
                'mapped': True,
                'chr': best_chain['ref_chr'],
                'pos': best_chain['ref_pos'],
                'score': best_chain['score']
            })
        else:
            results.append({
                'read_id': read['id'],
                'mapped': False,
                'chr': None,
                'pos': None
            })
    
    # Compare with ground truth
    print("\n" + "="*80)
    print("ACCURACY ANALYSIS")
    print("="*80)
    
    correct = 0
    wrong_pos = 0
    unmapped = 0
    tolerance = 10000  # 10kb
    
    for result in results:
        read_id = result['read_id']
        truth = ground_truth[read_id]
        
        if not result['mapped']:
            unmapped += 1
        elif result['chr'] == truth['chr']:
            if abs(result['pos'] - truth['start']) <= tolerance:
                correct += 1
            else:
                wrong_pos += 1
        else:
            # Wrong chromosome (should not happen with Chr22-only index!)
            wrong_pos += 1
    
    total = len(results)
    print(f"\nTotal reads: {total}")
    print(f"Correct (within 10kb): {correct} ({correct/total*100:.1f}%)")
    print(f"Wrong position: {wrong_pos} ({wrong_pos/total*100:.1f}%)")
    print(f"Unmapped: {unmapped} ({unmapped/total*100:.1f}%)")
    print()
    
    print("="*80)
    print("COMPARISON WITH PREVIOUS RESULTS")
    print("="*80)
    print()
    print("Full Genome Model:           Chr22 Focused Model:")
    print(f"  Correct: 43%               Correct: {correct/total*100:.1f}%")
    print(f"  Wrong: 49.5%               Wrong: {wrong_pos/total*100:.1f}%")
    print(f"  Unmapped: 7.5%             Unmapped: {unmapped/total*100:.1f}%")
    print()
    
    if correct / total >= 0.90:
        print("✅ SUCCESS: Accuracy >90%! Training helped significantly!")
    elif correct / total >= 0.70:
        print("⚠️  MODERATE: Accuracy improved but not to target")
    else:
        print("❌ FAILED: Accuracy still poor despite focused training")
    
    print()

if __name__ == '__main__':
    main()
