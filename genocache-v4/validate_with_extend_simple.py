#!/usr/bin/env python3
"""
Simple EXTEND Phase Validation

Uses working components + EXTEND phase to generate SAM output
Bypasses model complexity by using direct alignment approach
"""

import sys
import time
from pathlib import Path
import parasail

def load_fasta(fasta_file):
    """Load FASTA file"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def write_sam_header(sam_file):
    """Write minimal SAM header"""
    with open(sam_file, 'w') as f:
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        # Add main chromosomes
        chrs = [
            ("NC_000022.11", 50818468),  # chr22
            ("NC_000016.10", 90338345),  # chr16  
            ("NC_000013.11", 114364328), # chr13
            ("NC_000014.9", 107043718),  # chr14
        ]
        for chr_name, length in chrs:
            f.write(f"@SQ\tSN:{chr_name}\tLN:{length}\n")
        f.write("@PG\tID:genocache\tPN:genocache-extend\tVN:4.0\n")

def align_read_to_region(read_seq, ref_seq, ref_chr, ref_start, debug=False):
    """Align read to reference region using parasail"""
    matrix = parasail.matrix_create("ACGT", 2, -4)
    
    result = parasail.sg_qx_trace_scan_16(
        read_seq, ref_seq, 4, 2, matrix
    )
    
    if debug:
        print(f"      Debug: score={result.score}")
    
    # Lower threshold - any positive score
    if result.score < 10:  # Very low threshold for testing
        return None
    
    # Parse CIGAR
    cigar_str = result.cigar.decode if hasattr(result.cigar, 'decode') else str(len(read_seq)) + "M"
    
    # Try to get position
    ref_pos = getattr(result, 'ref_begin', 0)
    if ref_pos == 0 and hasattr(result, 'end_ref'):
        ref_pos = max(0, result.end_ref - len(read_seq))
    
    return {
        'chr': ref_chr,
        'pos': ref_start + ref_pos + 1,  # 1-based SAM
        'score': result.score,
        'cigar': cigar_str if cigar_str else f"{len(read_seq)}M",
        'mapq': 60 if result.score > 1000 else (40 if result.score > 500 else 20)
    }

def extend_phase_simple(read_id, read_seq, ref_sequences, candidates):
    """
    EXTEND phase: Align to each candidate, pick best by score
    
    This is THE FIX!
    """
    print(f"  [EXTEND] Aligning to {len(candidates)} candidates...")
    
    alignments = []
    
    for i, (chr_name, start, end) in enumerate(candidates):
        if chr_name not in ref_sequences:
            print(f"    {i+1}. {chr_name}: not loaded, skipping")
            continue
        
        # Extract reference region (expand search to ±5kb)
        ref_seq = ref_sequences[chr_name]
        region_start = max(0, start - 5000)
        region_end = min(len(ref_seq), end + 5000)
        ref_region = ref_seq[region_start:region_end]
        
        # Align
        alignment = align_read_to_region(read_seq, ref_region, chr_name, region_start, debug=True)
        
        if alignment:
            alignments.append(alignment)
            print(f"    {i+1}. {chr_name}:{start}-{end} → score={alignment['score']}, pos={alignment['pos']}")
        else:
            print(f"    {i+1}. {chr_name}:{start}-{end} → no alignment (score < 10)")
    
    if not alignments:
        print(f"  [EXTEND] No valid alignments found")
        return None
    
    # Pick best by ALIGNMENT SCORE (not seed count!)
    best = max(alignments, key=lambda x: x['score'])
    
    print(f"  [EXTEND] Best: {best['chr']}:{best['pos']} (score={best['score']})")
    
    return best

def write_sam_alignment(sam_file, read_id, read_seq, alignment):
    """Write alignment to SAM"""
    with open(sam_file, 'a') as f:
        if alignment is None:
            # Unmapped
            f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
        else:
            # Mapped
            flag = 0
            chr_name = alignment['chr']
            pos = alignment['pos']
            mapq = alignment['mapq']
            cigar = alignment['cigar']
            
            tags = f"AS:i:{alignment['score']}"
            
            f.write(f"{read_id}\t{flag}\t{chr_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read_seq}\t*\t{tags}\n")

def main():
    print("=" * 80)
    print("SIMPLE EXTEND PHASE VALIDATION")
    print("=" * 80)
    print()
    
    # Candidate chromosomes to test (based on OLD method failures)
    test_candidates_per_read = {
        'read_0': [
            ('NC_000016.10', 67964985, 67965985),  # OLD picked this (WRONG)
            ('NC_000022.11', 41905001, 41906001),  # Correct
        ],
        'read_5': [
            ('NC_000013.11', 99920954, 99921954),  # OLD picked this (WRONG)
            ('NC_000022.11', 23173001, 23174001),  # Correct
        ],
        'read_7': [
            ('NC_000014.9', 79197104, 79198104),    # OLD picked this (WRONG)
            ('NC_000022.11', 47701001, 47702001),   # Correct
        ],
        'read_9': [
            ('NC_000013.11', 46939113, 46940113),   # OLD picked this (WRONG)
            ('NC_000022.11', 39607001, 39608001),   # Correct
        ],
    }
    
    # Load reference
    print("Loading reference genome (main chromosomes)...")
    ref_file = '../GRCh38.fa'
    ref_sequences = load_fasta(ref_file)
    print(f"  Loaded {len(ref_sequences)} sequences")
    
    # Filter to test chromosomes
    test_chrs = ['NC_000022.11', 'NC_000016.10', 'NC_000013.11', 'NC_000014.9']
    ref_sequences = {k: v for k, v in ref_sequences.items() if k in test_chrs}
    print(f"  Using {len(ref_sequences)} test chromosomes")
    for chr_name, seq in ref_sequences.items():
        print(f"    {chr_name}: {len(seq):,} bp")
    print()
    
    # Load test reads
    print("Loading test reads...")
    test_reads_file = 'test_10_reads_exact.fa'
    all_reads = load_fasta(test_reads_file)
    test_reads = {k: v for k, v in all_reads.items() if k in test_candidates_per_read}
    print(f"  Loaded {len(test_reads)} test reads")
    print()
    
    # Initialize SAM output
    output_sam = 'test_with_extend_simple.sam'
    print(f"Initializing SAM output: {output_sam}")
    write_sam_header(output_sam)
    print(f"  ✅ SAM header written")
    print()
    
    print("=" * 80)
    print("RUNNING EXTEND PHASE")
    print("=" * 80)
    print()
    
    total_time = 0
    results = []
    
    for read_id in sorted(test_reads.keys()):
        read_seq = test_reads[read_id]
        candidates = test_candidates_per_read[read_id]
        
        print(f"\n{read_id} ({len(read_seq)} bp)")
        print(f"  Candidates to test: {len(candidates)}")
        for i, (chr_name, start, end) in enumerate(candidates):
            print(f"    {i+1}. {chr_name}:{start}-{end}")
        
        start_time = time.time()
        
        # Run EXTEND phase
        alignment = extend_phase_simple(read_id, read_seq, ref_sequences, candidates)
        
        # Write to SAM
        write_sam_alignment(output_sam, read_id, read_seq, alignment)
        
        elapsed = time.time() - start_time
        total_time += elapsed
        
        results.append({
            'read_id': read_id,
            'mapped': alignment is not None,
            'chr': alignment['chr'] if alignment else None
        })
        
        print(f"  Time: {elapsed:.3f}s")
    
    print()
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    print()
    
    mapped = sum(1 for r in results if r['mapped'])
    print(f"Total reads: {len(results)}")
    print(f"Mapped: {mapped}/{len(results)} ({100*mapped/len(results):.1f}%)")
    print(f"Total time: {total_time:.3f}s")
    print(f"Average: {total_time/len(results):.3f}s per read")
    print()
    
    print(f"✅ SAM output written to: {output_sam}")
    print()
    print("Next: Compare with minimap2:")
    print(f"  python3 compare_chromosome_accuracy.py {output_sam} minimap2_same_10reads.sam")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
