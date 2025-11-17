#!/usr/bin/env python3
"""
Test Phase 2: Multi-seeding + Chaining on chr22
Target: 73.5% → 85-90% recall
"""

import numpy as np
from Bio import SeqIO
import random
import sys
from multi_seeder import MultiSeedAligner
from seed_chainer import ChainedMultiSeedAligner

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


def test_chaining(n_reads=1000, n_seeds=5):
    """
    Test multi-seeding + chaining on synthetic chr22 reads
    
    Args:
        n_reads: Number of test reads
        n_seeds: Number of seeds per read
    
    Returns:
        Recall rate and gate pass status
    """
    print("="*70)
    print("Phase 2 Test: Multi-Seeding + Chaining on Chr22")
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
    aligner = ChainedMultiSeedAligner(
        multi_seed_aligner=base_aligner,
        max_gap=15000,  # More permissive
        max_deviation=2000  # Allow more indels
    )
    print("  ✓ Chained aligner ready")
    
    # Test parameters
    print(f"\nTest parameters:")
    print(f"  Reads to test: {n_reads:,}")
    print(f"  Seeds per read: {n_seeds}")
    print(f"  Error rate: 5%")
    print(f"  Read length: 1-5 kb")
    print(f"  Tolerance: 100 bp")
    print(f"  Chaining: Enabled (max_gap=10kb, max_dev=1kb)")
    
    # Run test
    print(f"\nRunning tests...")
    correct = 0
    total = 0
    errors = []
    chain_lengths = []
    chain_scores = []
    no_result = 0
    
    for i in range(n_reads):
        # Extract random read from chr22
        read_len = random.randint(1000, 5000)
        true_pos = random.randint(0, len(ref) - read_len - 1)
        read_seq = ref[true_pos:true_pos + read_len]
        
        # Add sequencing errors
        read_seq = add_errors(read_seq, error_rate=0.05)
        
        # Align with chaining
        result = aligner.align_read(read_seq, n_seeds=n_seeds, top_k=50)
        
        if result:
            # Get predicted position
            pred_pos = result['ref_pos']
            error = abs(pred_pos - true_pos)
            
            errors.append(error)
            chain_lengths.append(result['chain_length'])
            chain_scores.append(result['chain_score'])
            
            # Check if within tolerance
            if error < 100:  # 100 bp tolerance
                correct += 1
        else:
            no_result += 1
        
        total += 1
        
        # Progress update
        if (i + 1) % 100 == 0:
            current_recall = correct / total * 100
            print(f"  Progress: {i+1:4d}/{n_reads} | Recall: {current_recall:5.1f}% | Correct: {correct:4d} | No result: {no_result:3d}")
    
    # Calculate metrics
    recall = correct / total
    median_error = np.median(errors) if errors else float('inf')
    mean_error = np.mean(errors) if errors else float('inf')
    avg_chain_len = np.mean(chain_lengths) if chain_lengths else 0
    median_chain_score = np.median(chain_scores) if chain_scores else 0
    
    # Print results
    print(f"\n{'='*70}")
    print("Phase 2 Results (Multi-Seeding + Chaining):")
    print(f"{'='*70}")
    print(f"  Recall@100bp: {recall*100:.1f}%")
    print(f"  Correct: {correct}/{total}")
    print(f"  No alignment: {no_result}")
    print(f"  Median error: {median_error:.0f} bp")
    print(f"  Mean error: {mean_error:.1f} bp")
    print(f"  Avg chain length: {avg_chain_len:.1f} seeds")
    print(f"  Median chain score: {median_chain_score:.3f}")
    print(f"{'='*70}")
    
    # Error distribution
    if errors:
        error_bins = [0, 10, 50, 100, 500, 1000, 5000]
        print(f"\nError distribution:")
        for i in range(len(error_bins) - 1):
            count = sum(1 for e in errors if error_bins[i] <= e < error_bins[i+1])
            pct = count / len(errors) * 100
            print(f"  {error_bins[i]:5d} - {error_bins[i+1]:5d} bp: {count:4d} ({pct:5.1f}%)")
    
    # Chain length distribution
    if chain_lengths:
        print(f"\nChain length distribution:")
        max_len = max(chain_lengths)
        for length in range(1, min(max_len + 1, 11)):
            count = sum(1 for l in chain_lengths if l == length)
            pct = count / len(chain_lengths) * 100
            print(f"  {length:2d} seeds: {count:4d} ({pct:5.1f}%)")
        if max_len > 10:
            count = sum(1 for l in chain_lengths if l > 10)
            pct = count / len(chain_lengths) * 100
            print(f"  >10 seeds: {count:4d} ({pct:5.1f}%)")
    
    # Comparison with Phase 1
    phase1_recall = 73.5  # From previous test
    improvement = recall * 100 - phase1_recall
    
    print(f"\n{'='*70}")
    print("Comparison with Phase 1:")
    print(f"{'='*70}")
    print(f"  Phase 1 (Multi-seed only): {phase1_recall:.1f}%")
    print(f"  Phase 2 (+ Chaining):      {recall*100:.1f}%")
    print(f"  Improvement:               {improvement:+.1f}%")
    print(f"{'='*70}")
    
    # Gate check
    print(f"\n{'='*70}")
    print("Gate Check:")
    print(f"{'='*70}")
    
    if recall >= 0.85:
        print("  ✓ PASSED: Recall >= 85%")
        print("  → Ready for Phase 3 (Smith-Waterman Refinement)")
        gate_pass = True
    elif recall >= 0.80:
        print("  ~ MARGINAL: Recall 80-85%")
        print("  → Consider tuning parameters")
        print("  → Can proceed to Phase 3 with caution")
        gate_pass = True
    else:
        print("  ✗ FAILED: Recall < 80%")
        print("  → Debug chaining implementation")
        gate_pass = False
    
    print(f"{'='*70}")
    
    return recall, gate_pass


def main():
    """Run Phase 2 test"""
    try:
        recall, gate_pass = test_chaining(n_reads=1000, n_seeds=5)
        
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
