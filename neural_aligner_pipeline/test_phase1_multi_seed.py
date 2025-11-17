#!/usr/bin/env python3
"""
Test Phase 1: Multi-seeding improvement on chr22
Target: 70.6% → 80-85% recall
"""

import numpy as np
from Bio import SeqIO
import random
import sys
from multi_seeder import MultiSeedAligner

def add_errors(seq, error_rate=0.05):
    """Add substitution errors to simulate sequencing errors"""
    seq = list(seq)
    bases = ['A', 'C', 'G', 'T']
    
    for i in range(len(seq)):
        if random.random() < error_rate:
            # Substitute with different base
            seq[i] = random.choice([b for b in bases if b != seq[i]])
    
    return ''.join(seq)


def load_chr22_reference():
    """Load chr22 from reference genome"""
    print("Loading chr22 reference...")
    
    for record in SeqIO.parse("GRCh38.fa", "fasta"):
        # Chr22 can be named different ways
        if "NC_000022" in record.id or record.id == "chr22":
            ref = str(record.seq).upper()
            print(f"  ✓ Chr22 loaded: {len(ref):,} bp")
            return ref
    
    print("  ✗ Chr22 not found in reference!")
    return None


def test_multi_seed(n_reads=1000, n_seeds=5):
    """
    Test multi-seeding on synthetic chr22 reads
    
    Args:
        n_reads: Number of test reads
        n_seeds: Number of seeds per read
    
    Returns:
        Recall rate
    """
    print("="*70)
    print("Phase 1 Test: Multi-Seeding on Chr22")
    print("="*70)
    
    # Load reference
    ref = load_chr22_reference()
    if not ref:
        return 0.0
    
    # Create aligner
    print("\nInitializing multi-seed aligner...")
    aligner = MultiSeedAligner(
        encoder_path="nal_encoder_best.pt",
        index_path="faiss_index_improved.idx",
        positions_path="ref_positions_improved.npy"
    )
    
    # Test parameters
    print(f"\nTest parameters:")
    print(f"  Reads to test: {n_reads:,}")
    print(f"  Seeds per read: {n_seeds}")
    print(f"  Error rate: 5%")
    print(f"  Read length: 1-5 kb")
    print(f"  Tolerance: 100 bp")
    
    # Run test
    print(f"\nRunning tests...")
    correct = 0
    total = 0
    errors = []
    scores = []
    n_seeds_matched = []
    
    for i in range(n_reads):
        # Extract random read from chr22
        read_len = random.randint(1000, 5000)
        true_pos = random.randint(0, len(ref) - read_len - 1)
        read_seq = ref[true_pos:true_pos + read_len]
        
        # Add sequencing errors
        read_seq = add_errors(read_seq, error_rate=0.05)
        
        # Align with multi-seeding
        results = aligner.align_read(read_seq, n_seeds=n_seeds, top_k=20)
        
        if results:
            # Take top result
            pred_pos = results[0]['ref_pos']
            error = abs(pred_pos - true_pos)
            
            errors.append(error)
            scores.append(results[0]['total_score'])
            n_seeds_matched.append(results[0]['n_seeds'])
            
            # Check if within tolerance
            if error < 100:  # 100 bp tolerance
                correct += 1
        
        total += 1
        
        # Progress update
        if (i + 1) % 100 == 0:
            current_recall = correct / total * 100
            print(f"  Progress: {i+1:4d}/{n_reads} | Recall: {current_recall:5.1f}% | Correct: {correct:4d}")
    
    # Calculate metrics
    recall = correct / total
    median_error = np.median(errors) if errors else float('inf')
    mean_error = np.mean(errors) if errors else float('inf')
    median_score = np.median(scores) if scores else 0
    avg_seeds = np.mean(n_seeds_matched) if n_seeds_matched else 0
    
    # Print results
    print(f"\n{'='*70}")
    print("Phase 1 Results (Multi-Seeding):")
    print(f"{'='*70}")
    print(f"  Recall@100bp: {recall*100:.1f}%")
    print(f"  Correct: {correct}/{total}")
    print(f"  Median error: {median_error:.0f} bp")
    print(f"  Mean error: {mean_error:.1f} bp")
    print(f"  Median score: {median_score:.3f}")
    print(f"  Avg seeds matched: {avg_seeds:.1f}/{n_seeds}")
    print(f"{'='*70}")
    
    # Error distribution
    if errors:
        error_bins = [0, 10, 50, 100, 500, 1000, 5000]
        print(f"\nError distribution:")
        for i in range(len(error_bins) - 1):
            count = sum(1 for e in errors if error_bins[i] <= e < error_bins[i+1])
            pct = count / len(errors) * 100
            print(f"  {error_bins[i]:5d} - {error_bins[i+1]:5d} bp: {count:4d} ({pct:5.1f}%)")
    
    # Gate check
    print(f"\n{'='*70}")
    print("Gate Check:")
    print(f"{'='*70}")
    
    if recall >= 0.80:
        print("  ✓ PASSED: Recall >= 80%")
        print("  → Ready for Phase 2 (Chaining)")
        gate_pass = True
    elif recall >= 0.75:
        print("  ~ MARGINAL: Recall 75-80%")
        print("  → Tune parameters before Phase 2")
        gate_pass = False
    else:
        print("  ✗ FAILED: Recall < 75%")
        print("  → Debug multi-seeding implementation")
        gate_pass = False
    
    print(f"{'='*70}")
    
    return recall, gate_pass


def main():
    """Run Phase 1 test"""
    try:
        recall, gate_pass = test_multi_seed(n_reads=1000, n_seeds=5)
        
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
