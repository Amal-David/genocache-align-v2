#!/usr/bin/env python3
"""
GenoCache Production Aligner - FASTQ Input

Align FASTQ reads and output SAM format
Ready for GIAB validation
"""

import argparse
import numpy as np
import torch
import faiss
import edlib
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict
import time

from genocache_encoder import GenoCacheEncoder


def seq_to_tokens(seq):
    """Convert DNA sequence to integer tokens"""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    return [mapping.get(ch.upper(), 4) for ch in seq]


def load_fasta(fasta_path):
    """Load reference sequences"""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def load_positions(positions_file):
    """Load position mappings"""
    positions = []
    with open(positions_file, 'r') as f:
        for line in f:
            chrom, pos = line.strip().split('\t')
            positions.append((chrom, int(pos)))
    return positions


def sparse_dp_chain(hits, bandwidth=10000, gap_penalty=0.01, max_skip=100):
    """Sparse DP chaining"""
    if not hits:
        return [], 0.0
    
    n = len(hits)
    hits = sorted(hits, key=lambda x: x[0])
    
    dp = np.zeros(n, dtype=np.float32)
    parent = np.full(n, -1, dtype=np.int32)
    
    for i in range(n):
        dp[i] = hits[i][2]
    
    for i in range(n):
        query_i, ref_i, sim_i = hits[i]
        start_j = max(0, i - max_skip)
        
        for j in range(start_j, i):
            query_j, ref_j, sim_j = hits[j]
            
            query_gap = query_i - query_j
            ref_gap = ref_i - ref_j
            
            if query_gap <= 0 or ref_gap <= 0:
                continue
            
            if ref_gap > bandwidth:
                continue
            
            gap_diff = abs(ref_gap - query_gap)
            gap_cost = gap_penalty * gap_diff
            
            score = dp[j] + sim_i - gap_cost
            
            if score > dp[i]:
                dp[i] = score
                parent[i] = j
    
    best_i = np.argmax(dp)
    
    chain = []
    i = best_i
    while i != -1:
        chain.append(hits[i])
        i = parent[i]
    
    chain.reverse()
    return chain, float(dp[best_i])


def find_best_chain(seed_matches, bandwidth=10000, gap_penalty=0.01, max_skip=100):
    """Find best chain across all chromosomes"""
    by_chrom = defaultdict(list)
    for seed_off, chrom, ref_pos, sim in seed_matches:
        by_chrom[chrom].append((seed_off, ref_pos, sim))
    
    chains = []
    for chrom, hits in by_chrom.items():
        if len(hits) < 1:
            continue
        
        chain, score = sparse_dp_chain(hits, bandwidth, gap_penalty, max_skip)
        
        if chain:
            chains.append({
                'chrom': chrom,
                'chain': chain,
                'score': score
            })
    
    if not chains:
        return None
    
    chains.sort(key=lambda x: -x['score'])
    
    result = chains[0]
    result['second_best_score'] = chains[1]['score'] if len(chains) > 1 else 0.0
    
    return result


def extend_alignment_with_edlib(query_seq, ref_seq, seed_chain, ref_start_pos):
    """Extend seed chain to full read alignment using edlib"""
    if not seed_chain:
        return None
    
    query_len = len(query_seq)
    first_seed_offset = seed_chain[0][0]
    first_ref_pos = seed_chain[0][1]
    
    read_start = first_ref_pos - first_seed_offset
    read_start = max(0, read_start)
    
    padding = 500
    ref_region_start = max(0, read_start - padding)
    ref_region_end = min(len(ref_seq), read_start + query_len + padding)
    
    ref_region = ref_seq[ref_region_start:ref_region_end]
    
    result = edlib.align(
        query_seq,
        ref_region,
        mode="HW",
        task="path",
        k=-1
    )
    
    if result['editDistance'] == -1:
        return {
            'ref_start': read_start,
            'ref_end': read_start + query_len,
            'cigar': f"{query_len}M",
            'alignment_score': 0,
            'edit_distance': -1
        }
    
    cigar_tuples = result.get('cigar', None)
    if cigar_tuples:
        cigar = edlib_to_sam_cigar(cigar_tuples)
    else:
        cigar = f"{query_len}M"
    
    actual_ref_start = ref_region_start + result['locations'][0][0]
    actual_ref_end = ref_region_start + result['locations'][0][1] + 1
    
    return {
        'ref_start': actual_ref_start,
        'ref_end': actual_ref_end,
        'cigar': cigar,
        'alignment_score': len(query_seq) - result['editDistance'],
        'edit_distance': result['editDistance']
    }


def edlib_to_sam_cigar(edlib_cigar):
    """Convert edlib CIGAR to SAM CIGAR format"""
    if not edlib_cigar:
        return ""
    
    sam_cigar = []
    i = 0
    while i < len(edlib_cigar):
        if edlib_cigar[i].isdigit():
            count_str = ""
            while i < len(edlib_cigar) and edlib_cigar[i].isdigit():
                count_str += edlib_cigar[i]
                i += 1
            count = int(count_str)
            
            if i < len(edlib_cigar):
                op = edlib_cigar[i]
                i += 1
                
                if op == '=' or op == 'X':
                    sam_op = 'M'
                elif op == 'I':
                    sam_op = 'I'
                elif op == 'D':
                    sam_op = 'D'
                else:
                    sam_op = 'M'
                
                sam_cigar.append(f"{count}{sam_op}")
        else:
            i += 1
    
    return ''.join(sam_cigar) if sam_cigar else "*"


def calculate_mapq(best_score, second_best_score, chain_length):
    """Calculate MAPQ"""
    if chain_length == 0:
        return 0
    
    score_diff = best_score - second_best_score
    
    if score_diff <= 0:
        return 0
    
    uniqueness = min(40, int(20 * score_diff))
    length_bonus = min(20, int(5 * np.log1p(chain_length)))
    
    mapq = uniqueness + length_bonus
    mapq = max(0, min(60, mapq))
    
    return mapq


def align_read(
    model, index, positions, ref_sequences,
    read_id, read_sequence,
    seed_len=512,
    num_seeds=15,
    k_per_seed=50,
    bandwidth=10000,
    gap_penalty=0.01,
    max_skip=100,
    device='cuda'
):
    """Align a single read"""
    if len(read_sequence) < seed_len:
        return None
    
    # Extract seeds
    if len(read_sequence) == seed_len:
        offsets = [0]
    else:
        stride = (len(read_sequence) - seed_len) // (num_seeds - 1)
        offsets = [i * stride for i in range(num_seeds)]
    
    seed_matches = []
    
    for offset in offsets:
        seed = read_sequence[offset:offset + seed_len]
        
        if seed.count('N') > seed_len * 0.1:
            continue
        
        tokens = seq_to_tokens(seed)
        tokens_tensor = torch.tensor([tokens], dtype=torch.long, device=device)
        
        with torch.no_grad():
            emb = model(tokens_tensor).cpu().numpy()
        
        faiss.normalize_L2(emb)
        D, I = index.search(emb, k_per_seed)
        
        for idx, similarity in zip(I[0], D[0]):
            chrom, pos = positions[idx]
            seed_matches.append((offset, chrom, pos, similarity))
    
    if not seed_matches:
        return None
    
    # Find best chain
    best_chain_info = find_best_chain(seed_matches, bandwidth, gap_penalty, max_skip)
    
    if not best_chain_info:
        return None
    
    chrom = best_chain_info['chrom']
    seed_chain = best_chain_info['chain']
    chain_score = best_chain_info['score']
    second_best_score = best_chain_info['second_best_score']
    
    if chrom not in ref_sequences:
        return None
    
    ref_seq = ref_sequences[chrom]
    
    # Extend alignment
    first_seed_ref_pos = seed_chain[0][1]
    alignment = extend_alignment_with_edlib(
        read_sequence,
        ref_seq,
        seed_chain,
        first_seed_ref_pos
    )
    
    if not alignment:
        return None
    
    # Calculate MAPQ
    mapq = calculate_mapq(chain_score, second_best_score, len(seed_chain))
    
    return {
        'read_id': read_id,
        'chrom': chrom,
        'pos': alignment['ref_start'] + 1,  # 1-based
        'mapq': mapq,
        'cigar': alignment['cigar'],
        'sequence': read_sequence,
        'aligned': True,
        'chain_length': len(seed_chain),
        'edit_distance': alignment['edit_distance'],
        'alignment_score': alignment['alignment_score']
    }


def write_sam_header(f, ref_sequences):
    """Write SAM header"""
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    for chrom, seq in ref_sequences.items():
        f.write(f"@SQ\tSN:{chrom}\tLN:{len(seq)}\n")
    f.write("@PG\tID:genocache\tPN:GenoCache\tVN:3.0\tCL:genocache_fastq\n")


def write_sam_record(f, alignment):
    """Write SAM alignment record"""
    if not alignment or not alignment.get('aligned', False):
        f.write(f"{alignment['read_id']}\t4\t*\t0\t0\t*\t*\t0\t0\t{alignment['sequence']}\t*\n")
        return
    
    flag = 0
    
    f.write(f"{alignment['read_id']}\t{flag}\t{alignment['chrom']}\t{alignment['pos']}\t"
           f"{alignment['mapq']}\t{alignment['cigar']}\t*\t0\t0\t"
           f"{alignment['sequence']}\t*\t"
           f"AS:i:{alignment['alignment_score']}\t"
           f"NM:i:{alignment['edit_distance']}\t"
           f"cs:i:{alignment['chain_length']}\n")


def main():
    parser = argparse.ArgumentParser(description="GenoCache FASTQ Aligner")
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output-sam", required=True, help="Output SAM file")
    parser.add_argument("--num-seeds", type=int, default=15)
    parser.add_argument("--k-per-seed", type=int, default=50)
    parser.add_argument("--bandwidth", type=int, default=10000)
    parser.add_argument("--gap-penalty", type=float, default=0.01)
    parser.add_argument("--max-skip", type=int, default=100)
    
    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print("="*80)
    print("GenoCache Production Aligner")
    print("="*80)
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
    
    # Load positions
    print("Loading positions...")
    positions = load_positions(args.positions)
    print(f"  ✅ Positions: {len(positions):,}")
    
    # Load reference
    print("Loading reference...")
    ref_sequences = load_fasta(args.fasta)
    print(f"  ✅ Reference: {len(ref_sequences)} sequences")
    print()
    
    # Load reads
    print(f"Loading reads from {args.fastq}...")
    reads = []
    for record in SeqIO.parse(args.fastq, "fastq"):
        reads.append((record.id, str(record.seq)))
    print(f"  ✅ Loaded {len(reads)} reads")
    print()
    
    # Align
    print("Aligning reads...")
    start_time = time.time()
    
    with open(args.output_sam, 'w') as sam_out:
        # Write header
        write_sam_header(sam_out, ref_sequences)
        
        # Align reads
        for read_id, read_seq in tqdm(reads, desc="Reads"):
            alignment = align_read(
                model, index, positions, ref_sequences,
                read_id, read_seq,
                num_seeds=args.num_seeds,
                k_per_seed=args.k_per_seed,
                bandwidth=args.bandwidth,
                gap_penalty=args.gap_penalty,
                max_skip=args.max_skip,
                device=device
            )
            
            if alignment is None:
                alignment = {'read_id': read_id, 'sequence': read_seq, 'aligned': False}
            
            write_sam_record(sam_out, alignment)
    
    elapsed = time.time() - start_time
    
    print()
    print(f"✅ Alignment complete!")
    print(f"   Time: {elapsed:.1f}s")
    print(f"   Speed: {len(reads)/elapsed:.1f} reads/sec")
    print(f"   Output: {args.output_sam}")
    print()


if __name__ == "__main__":
    main()
