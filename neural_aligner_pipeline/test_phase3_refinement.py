#!/usr/bin/env python3
"""
Test Phase 3: Multi-seeding + Smith-Waterman Refinement
Expected improvement: 73.5% → 90-95%
"""

import numpy as np
from Bio import SeqIO
import random
import sys
from multi_seeder import MultiSeedAligner
from seed_chainer import ChainedMultiSeedAligner
from refine_alignment import AlignmentRefiner, RefinedAligner

def add_errors(seq, error_rate=0.05):
    """Add substitution errors to simulate sequencing errors"""
    seq = list(seq)
    bases = ['A', 'C', 'G', 'T']
    
    for i in range(len(seq)):
        if random.random() < error_rate:
            seq[i] = random.choice([b for b in bases if b != seq[i]])
    
    return ''.join(seq)


def load_chr22_reference():
    """Load chr22 from reference genome"""
    print("Loading chr22 reference...")
    
    for record in SeqIO.parse("GRCh38.fa", "fasta"):
        if "NC_000022" in record.id or record.id == "chr22":
            ref = str(record.seq).upper()
            print(f"  ✓ Chr22 loaded: {len(ref):,} bp")
            return ref
    
    print("  ✗ Chr22 not found in reference!")
    return None


def test_sw_refinement(n_reads=1000, n_seeds=5):
    """
    Test multi-seeding + SW refinement on synthetic chr22 reads
    
    Args:
        n_reads: Number of test reads
        n_seeds: Number of seeds per read
    
    Returns:
        Recall rate and gate pass status
    """
    print("="*70)
    print("Phase 3 Test: Multi-Seeding + Smith-Waterman Refinement")
    print("="*70)
    
    # Load reference
    ref = load_chr22_reference()
    if not ref:
        return 0.0, False
    
    # Create base aligner
    print("\nInitializing multi-seed aligner...")
    base_aligner = MultiSeedAligner(
        encoder_path="nal_encoder_best.pt",
        index_path="faiss_index_improved.idx",
        positions_path="ref_positions_improved.npy"
    )
    
    # Create chained aligner
    print("Initializing seed chainer...")
    chained_aligner = ChainedMultiSeedAligner(
        multi_seed_aligner=base_aligner,
        max_gap=15000,
        max_deviation=2000
    )
    
    # Create SW refiner
    print("Loading Smith-Waterman refiner...")
    refiner = AlignmentRefiner(ref_fa="GRCh38.fa")
    
    # Create complete refined aligner
    print("Creating refined aligner...")
    aligner = RefinedAligner(
        chained_aligner=chained_aligner,
        refiner=refiner
    )
    print("  ✓ Complete pipeline ready")
    
    # Test parameters
    print(f"\nTest parameters:")
    print(f"  Reads to test: {n_reads:,}")
    print(f"  Seeds per read: {n_seeds}")
    print(f"  Error rate: 5%")
    print(f"  Read length: 1-5 kb")
    print(f"  Tolerance: 100 bp")
    print(f"  Pipeline: Multi-seed → Chain → SW Refinement")
    
    # Run test
    print(f"\nRunning tests...")
    correct = 0
    total = 0
    errors = []
    sw_scores = []
    identities = []
    refined_count = 0
    no_result = 0
    
    for i in range(n_reads):
        # Extract random read from chr22
        read_len = random.randint(1000, 5000)
        
        # Avoid regions with too many Ns
        max_attempts = 10
        for attempt in range(max_attempts):
            true_pos = random.randint(0, len(ref) - read_len - 1)
            read_seq = ref[true_pos:true_pos + read_len]
            
            if read_seq.count('N') < read_len * 0.2:  # Less than 20% Ns
                break
        
        # Add sequencing errors
        read_seq = add_errors(read_seq, error_rate=0.05)
        
        # Align with SW refinement
        result = aligner.align_read(read_seq, n_seeds=n_seeds)
        
        if result:
            # Get predicted position
            pred_pos = result['position']
            error = abs(pred_pos - true_pos)
            
            errors.append(error)
            
            if result.get('refined', False):
                refined_count += 1
                if 'score' in result:
                    sw_scores.append(result['score'])
                if 'identity' in result:
                    identities.append(result['identity'])
            
            # Check if within tolerance
            if error < 100:  # 100 bp tolerance
                correct += 1
        else:
            no_result += 1
        
        total += 1
        
        # Progress update
        if (i + 1) % 100 == 0:
            current_recall = correct / total * 100
            refined_pct = refined_count / total * 100
            print(f"  Progress: {i+1:4d}/{n_reads} | Recall: {current_recall:5.1f}% | Refined: {refined_pct:4.1f}% | No result: {no_result:3d}")
    
    # Calculate metrics
    recall = correct / total
    median_error = np.median(errors) if errors else float('inf')
    mean_error = np.mean(errors) if errors else float('inf')
    median_score = np.median(sw_scores) if sw_scores else 0
    avg_identity = np.mean(identities) if identities else 0
    refined_pct = refined_count / total * 100
    
    # Print results
    print(f"\n{'='*70}")
    print("Phase 3 Results (Multi-Seed + SW Refinement):")
    print(f"{'='*70}")
    print(f"  Recall@100bp: {recall*100:.1f}%")
    print(f"  Correct: {correct}/{total}")
    print(f"  No alignment: {no_result}")
    print(f"  Refined with SW: {refined_count}/{total} ({refined_pct:.1f}%)")
    print(f"  Median error: {median_error:.0f} bp")
    print(f"  Mean error: {mean_error:.1f} bp")
    print(f"  Median SW score: {median_score:.0f}")
    print(f"  Avg identity: {avg_identity*100:.1f}%")
    print(f"{'='*70}")
    
    # Error distribution
    if errors:
        error_bins = [0, 5, 10, 50, 100, 500, 1000]
        print(f"\nError distribution:")
        for i in range(len(error_bins) - 1):
            count = sum(1 for e in errors if error_bins[i] <= e < error_bins[i+1])
            pct = count / len(errors) * 100
            print(f"  {error_bins[i]:5d} - {error_bins[i+1]:5d} bp: {count:4d} ({pct:5.1f}%)")
        
        # Large errors
        large_errors = sum(1 for e in errors if e >= 1000)
        if large_errors > 0:
            pct = large_errors / len(errors) * 100
            print(f"  ≥1000 bp:           {large_errors:4d} ({pct:5.1f}%)")
    
    # Comparison with previous phases
    phase1_recall = 73.5
    phase2_recall = 73.5
    improvement_from_p1 = recall * 100 - phase1_recall
    improvement_from_p2 = recall * 100 - phase2_recall
    
    print(f"\n{'='*70}")
    print("Comparison with Previous Phases:")
    print(f"{'='*70}")
    print(f"  Phase 1 (Multi-seed):         {phase1_recall:.1f}%")
    print(f"  Phase 2 (+ Chaining):         {phase2_recall:.1f}%")
    print(f"  Phase 3 (+ SW Refinement):    {recall*100:.1f}%")
    print(f"  Improvement from Phase 1:     {improvement_from_p1:+.1f}%")
    print(f"  Improvement from Phase 2:     {improvement_from_p2:+.1f}%")
    print(f"{'='*70}")
    
    # Gate check
    print(f"\n{'='*70}")
    print("Gate Check:")
    print(f"{'='*70}")
    
    if recall >= 0.90:
        print("  ✓ PASSED: Recall >= 90%")
        print("  → Ready for Phase 4 (Quality Scores + BAM Output)")
        gate_pass = True
    elif recall >= 0.85:
        print("  ~ MARGINAL: Recall 85-90%")
        print("  → Consider tuning SW parameters")
        print("  → Can proceed to Phase 4 with caution")
        gate_pass = True
    elif recall >= 0.80:
        print("  ~ IMPROVEMENT: Recall 80-85%")
        print("  → Significant gain from SW")
        print("  → Proceed to Phase 4")
        gate_pass = True
    else:
        print("  ✗ NEEDS WORK: Recall < 80%")
        print("  → Debug SW integration")
        print("  → Check margin size and parameters")
        gate_pass = False
    
    print(f"{'='*70}")
    
    return recall, gate_pass


def main():
    """Run Phase 3 test"""
    try:
        recall, gate_pass = test_sw_refinement(n_reads=1000, n_seeds=5)
        
        # Exit code for scripting
        sys.exit(0 if gate_pass else 1)
        
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
