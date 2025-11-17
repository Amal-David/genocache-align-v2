#!/usr/bin/env python3
"""
Read alignment with SPARSE DP CHAINING (minimap2-style)

IMPROVEMENT OVER GREEDY CHAINING:
- Greedy: Extends from each seed, stops if not co-linear (chain_length = 1.1)
- Sparse DP: Finds optimal chain globally (expected chain_length = 8-15)

Algorithm:
1. Group hits by chromosome (same as before)
2. For each chromosome, use sparse DP:
   - dp[i] = best score ending at seed i
   - dp[i] = max(score[i], max(dp[j] + score[i] - gap_cost) for all valid j < i)
   - Valid = co-linear + within bandwidth + gap cost acceptable
3. Pick best chain across all chromosomes

Expected improvement: 42.6% → 75-85% accuracy
"""

import argparse
import json
import numpy as np
import torch
import faiss
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict

from genocache_encoder import GenoCacheEncoder


def seq_to_tokens(seq):
    """Convert DNA sequence to integer tokens"""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    return [mapping.get(ch.upper(), 4) for ch in seq]


def load_fasta(fasta_path):
    """Load reference sequences"""
    sequences = {}
    print(f"Loading reference: {fasta_path}")
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    print(f"  Loaded {len(sequences)} sequences")
    return sequences


def load_positions(positions_file):
    """Load position mappings"""
    positions = []
    with open(positions_file, 'r') as f:
        for line in f:
            chrom, pos = line.strip().split('\t')
            positions.append((chrom, int(pos)))
    return positions


def sparse_dp_chain(hits, bandwidth=5000, gap_penalty=0.01, max_skip=25, debug=False):
    """
    Sparse DP chaining (minimap2-style)
    
    Key improvements over greedy:
    - Considers ALL valid predecessors (not just next seed)
    - Bandwidth constraint: Only seeds within bandwidth on reference
    - Gap penalties: Penalize insertions/deletions
    - Max skip: Only look back max_skip seeds (sparse)
    
    Args:
        hits: [(seed_offset, ref_pos, similarity), ...] - all on same chromosome
        bandwidth: Maximum ref distance to consider (5kb)
        gap_penalty: Penalty per bp of gap difference (0.01)
        max_skip: Maximum number of seeds to look back (25)
        debug: Print debug info
        
    Returns:
        best_chain: [list of (seed_offset, ref_pos, similarity)]
    """
    if not hits:
        return []
    
    n = len(hits)
    
    # Sort by query position (seed offset)
    hits = sorted(hits, key=lambda x: x[0])
    
    # DP arrays
    dp = np.zeros(n, dtype=np.float32)  # Best score ending at i
    parent = np.full(n, -1, dtype=np.int32)  # Backtrack pointer
    
    # Initialize: each seed can start a chain
    for i in range(n):
        dp[i] = hits[i][2]  # similarity score
    
    # DP: Find best predecessor for each seed
    for i in range(n):
        query_i, ref_i, sim_i = hits[i]
        
        # Look back at most max_skip seeds
        start_j = max(0, i - max_skip)
        
        for j in range(start_j, i):
            query_j, ref_j, sim_j = hits[j]
            
            # Check co-linearity
            query_gap = query_i - query_j
            ref_gap = ref_i - ref_j
            
            # Both must be positive (monotonic)
            if query_gap <= 0 or ref_gap <= 0:
                continue
            
            # Bandwidth constraint: ref_gap must be reasonable
            if ref_gap > bandwidth:
                continue
            
            # Gap cost: penalize distance difference
            gap_diff = abs(ref_gap - query_gap)
            gap_cost = gap_penalty * gap_diff
            
            # Can we improve by extending from j?
            score = dp[j] + sim_i - gap_cost
            
            if score > dp[i]:
                dp[i] = score
                parent[i] = j
    
    # Find best endpoint
    best_i = np.argmax(dp)
    best_score = dp[best_i]
    
    # Backtrack to build chain
    chain = []
    i = best_i
    while i != -1:
        chain.append(hits[i])
        i = parent[i]
    
    chain.reverse()
    
    if debug:
        print(f"\n[SPARSE DP CHAINING]")
        print(f"  Input seeds: {n}")
        print(f"  Chain length: {len(chain)}")
        print(f"  Chain score: {best_score:.3f}")
        print(f"  Avg similarity: {np.mean([s for _, _, s in chain]):.3f}")
        if len(chain) >= 2:
            positions = [pos for _, pos, _ in chain]
            print(f"  Position span: {positions[0]:,} - {positions[-1]:,} ({positions[-1]-positions[0]:,} bp)")
    
    return chain


def find_best_chain_sparse_dp(seed_matches, bandwidth=5000, gap_penalty=0.01, max_skip=25, debug=False):
    """
    Find best chain across all chromosomes using sparse DP
    
    Args:
        seed_matches: [(seed_offset, chrom, ref_pos, similarity), ...]
        bandwidth: Maximum ref distance for co-linearity (5kb)
        gap_penalty: Penalty per bp of gap difference
        max_skip: Max seeds to look back
        debug: Print debug info
        
    Returns:
        best_chain: {
            'chrom': str,
            'start_pos': int,
            'chain_length': int,
            'chain_score': float,
            'seeds': [(seed_offset, ref_pos, similarity), ...]
        }
    """
    # Group by chromosome
    by_chrom = defaultdict(list)
    for seed_off, chrom, ref_pos, sim in seed_matches:
        by_chrom[chrom].append((seed_off, ref_pos, sim))
    
    if debug and not hasattr(find_best_chain_sparse_dp, '_debug_done'):
        find_best_chain_sparse_dp._debug_done = True
        print(f"\n[SPARSE DP - Per Chromosome]")
        print(f"  Total seed matches: {len(seed_matches)}")
        print(f"  Chromosomes: {len(by_chrom)}")
        print(f"  Top 5 chromosomes:")
        for chrom, hits in sorted(by_chrom.items(), key=lambda x: -len(x[1]))[:5]:
            print(f"    {chrom}: {len(hits)} hits")
    
    # Try chaining on each chromosome
    best_chain_info = None
    best_score = 0
    
    for chrom, hits in by_chrom.items():
        if len(hits) < 1:
            continue
        
        # Run sparse DP on this chromosome
        chain = sparse_dp_chain(
            hits,
            bandwidth=bandwidth,
            gap_penalty=gap_penalty,
            max_skip=max_skip,
            debug=False
        )
        
        if not chain:
            continue
        
        # Score chain
        chain_length = len(chain)
        avg_similarity = np.mean([s for _, _, s in chain])
        chain_score = chain_length * avg_similarity
        
        if chain_score > best_score:
            best_score = chain_score
            best_chain_info = {
                'chrom': chrom,
                'start_pos': chain[0][1],  # First ref position
                'chain_length': chain_length,
                'chain_score': chain_score,
                'avg_similarity': avg_similarity,
                'seeds': chain
            }
    
    if debug and best_chain_info:
        print(f"\n[BEST CHAIN]")
        print(f"  Chromosome: {best_chain_info['chrom']}")
        print(f"  Chain length: {best_chain_info['chain_length']} seeds")
        print(f"  Score: {best_chain_info['chain_score']:.3f}")
        print(f"  Avg similarity: {best_chain_info['avg_similarity']:.3f}")
    
    return best_chain_info


def align_read_with_sparse_dp(
    model, index, vectors_mmap, positions,
    read_sequence,
    seed_len=512,
    num_seeds=5,
    k_per_seed=50,
    bandwidth=5000,
    gap_penalty=0.01,
    max_skip=25,
    device='cuda'
):
    """
    Align read using sparse DP chaining
    """
    if len(read_sequence) < seed_len:
        return None
    
    # Calculate seed offsets
    if len(read_sequence) == seed_len:
        offsets = [0]
    else:
        stride = (len(read_sequence) - seed_len) // (num_seeds - 1)
        offsets = [i * stride for i in range(num_seeds)]
    
    # Extract and search seeds
    seed_matches = []
    
    for offset in offsets:
        seed = read_sequence[offset:offset + seed_len]
        
        if seed.count('N') > seed_len * 0.1:
            continue
        
        # Encode
        tokens = seq_to_tokens(seed)
        tokens_tensor = torch.tensor([tokens], dtype=torch.long, device=device)
        
        with torch.no_grad():
            emb = model(tokens_tensor).cpu().numpy()
        
        # Search
        faiss.normalize_L2(emb)
        D, I = index.search(emb, k_per_seed)
        
        # Collect matches
        for idx, similarity in zip(I[0], D[0]):
            chrom, pos = positions[idx]
            seed_matches.append((offset, chrom, pos, similarity))
    
    if not seed_matches:
        return None
    
    # Find best chain using sparse DP
    debug = not hasattr(align_read_with_sparse_dp, '_debug_done')
    if debug:
        align_read_with_sparse_dp._debug_done = True
    
    best_chain = find_best_chain_sparse_dp(
        seed_matches,
        bandwidth=bandwidth,
        gap_penalty=gap_penalty,
        max_skip=max_skip,
        debug=debug
    )
    
    if best_chain is None:
        return None
    
    # Calculate read start position
    first_seed_offset = best_chain['seeds'][0][0]
    first_ref_pos = best_chain['seeds'][0][1]
    read_start = first_ref_pos - first_seed_offset
    
    return {
        'chrom': best_chain['chrom'],
        'pos': read_start,
        'chain_length': best_chain['chain_length'],
        'chain_score': best_chain['chain_score'],
        'avg_similarity': best_chain['avg_similarity']
    }


def benchmark_alignment(args):
    """Benchmark alignment accuracy with sparse DP chaining"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print("="*80)
    print("GenoCache Alignment with SPARSE DP CHAINING")
    print("="*80)
    print()
    print("Algorithm improvements:")
    print("  - Sparse DP (considers all valid predecessors)")
    print(f"  - Bandwidth constraint: {args.bandwidth:,} bp")
    print(f"  - Gap penalty: {args.gap_penalty}")
    print(f"  - Max skip: {args.max_skip} seeds")
    print()
    print("Expected: 42.6% → 75-85% accuracy, 1.1 → 8-15 chain length")
    print()
    
    # Load model
    print("Loading model...")
    model = GenoCacheEncoder(seed_len=512, emb_dim=256)
    checkpoint = torch.load(args.checkpoint, map_location=device)
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model.load_state_dict(checkpoint)
    model = model.to(device)
    model.eval()
    print("  ✅ Model loaded")
    
    # Load index
    print("Loading index...")
    index = faiss.read_index(args.index)
    print(f"  ✅ Index: {index.ntotal:,} vectors")
    
    # Load vectors (not needed for IVF-Flat)
    vectors_mmap = None
    
    # Load positions
    print("Loading positions...")
    positions = load_positions(args.positions)
    print(f"  ✅ Positions: {len(positions):,}")
    
    # Load reference
    print("Loading reference...")
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    print(f"  ✅ Reference: {len(ref_seq):,} bp")
    print()
    
    # Generate reads
    print(f"Generating {args.num_reads} reads...")
    np.random.seed(42)
    read_len = args.read_length
    max_start = len(ref_seq) - read_len - 10
    
    reads = []
    for _ in range(args.num_reads * 2):
        start = np.random.randint(10, max_start)
        read_seq = ref_seq[start:start + read_len]
        if read_seq.count('N') < read_len * 0.05:
            reads.append((start, read_seq))
        if len(reads) >= args.num_reads:
            break
    
    print(f"  Generated {len(reads)} reads")
    print()
    
    # Align
    stats = {
        'total': 0,
        'mapped': 0,
        'correct': 0,
        'correct_chrom': 0
    }
    chain_lengths = []
    
    print("Aligning reads...")
    for true_pos, read_seq in tqdm(reads, desc="Reads"):
        stats['total'] += 1
        
        alignment = align_read_with_sparse_dp(
            model, index, vectors_mmap, positions,
            read_seq,
            num_seeds=args.num_seeds,
            k_per_seed=args.k_per_seed,
            bandwidth=args.bandwidth,
            gap_penalty=args.gap_penalty,
            max_skip=args.max_skip,
            device=device
        )
        
        if alignment is None:
            continue
        
        stats['mapped'] += 1
        chain_lengths.append(alignment['chain_length'])
        
        # Check correctness
        if alignment['chrom'] == args.chrom:
            stats['correct_chrom'] += 1
            if abs(alignment['pos'] - true_pos) <= args.tolerance:
                stats['correct'] += 1
    
    # Results
    print()
    print("="*80)
    print("RESULTS")
    print("="*80)
    print()
    print(f"Total reads: {stats['total']}")
    print(f"Mapped: {stats['mapped']} ({100*stats['mapped']/stats['total']:.1f}%)")
    print(f"Correct chromosome: {stats['correct_chrom']} ({100*stats['correct_chrom']/stats['total']:.1f}%)")
    print(f"Correctly mapped: {stats['correct']} ({100*stats['correct']/stats['total']:.1f}%)")
    print()
    print(f"⭐ ACCURACY: {100*stats['correct']/stats['total']:.2f}%")
    if stats['mapped'] > 0:
        print(f"   On-target rate: {100*stats['correct_chrom']/stats['mapped']:.1f}%")
    print()
    
    if chain_lengths:
        print("Chain length stats:")
        print(f"  Mean: {np.mean(chain_lengths):.1f} seeds")
        print(f"  Median: {np.median(chain_lengths):.0f} seeds")
        print(f"  Max: {np.max(chain_lengths)} seeds")
        print()
    
    # Compare with baseline
    print("="*80)
    print("IMPROVEMENT vs GREEDY BASELINE:")
    print("="*80)
    baseline_acc = 42.6
    baseline_chain = 1.1
    improvement_acc = (100*stats['correct']/stats['total']) - baseline_acc
    improvement_chain = np.mean(chain_lengths) - baseline_chain if chain_lengths else 0
    print(f"  Accuracy: {baseline_acc:.1f}% → {100*stats['correct']/stats['total']:.1f}% ({improvement_acc:+.1f}%)")
    if chain_lengths:
        print(f"  Chain length: {baseline_chain:.1f} → {np.mean(chain_lengths):.1f} seeds ({improvement_chain:+.1f})")
    print()
    print("="*80)
    
    # Save
    if args.output:
        results = {
            'stats': stats,
            'chain_lengths': [int(x) for x in chain_lengths],
            'baseline_comparison': {
                'baseline_accuracy': baseline_acc,
                'new_accuracy': 100*stats['correct']/stats['total'],
                'improvement': improvement_acc
            }
        }
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Alignment with sparse DP chaining")
    parser.add_argument("--checkpoint", required=True, help="Model checkpoint")
    parser.add_argument("--index", required=True, help="FAISS index")
    parser.add_argument("--positions", required=True, help="Position mappings")
    parser.add_argument("--fasta", required=True, help="Reference FASTA")
    parser.add_argument("--chrom", required=True, help="Chromosome to test")
    parser.add_argument("--num-reads", type=int, default=1000, help="Number of test reads")
    parser.add_argument("--read-length", type=int, default=2000, help="Read length")
    parser.add_argument("--num-seeds", type=int, default=5, help="Seeds per read")
    parser.add_argument("--k-per-seed", type=int, default=50, help="FAISS k per seed")
    parser.add_argument("--bandwidth", type=int, default=5000, help="Max ref gap for co-linearity (bp)")
    parser.add_argument("--gap-penalty", type=float, default=0.01, help="Penalty per bp of gap difference")
    parser.add_argument("--max-skip", type=int, default=25, help="Max seeds to look back in DP")
    parser.add_argument("--tolerance", type=int, default=1024, help="Correctness tolerance (bp)")
    parser.add_argument("--output", help="Output JSON file")
    
    args = parser.parse_args()
    benchmark_alignment(args)


if __name__ == "__main__":
    main()
