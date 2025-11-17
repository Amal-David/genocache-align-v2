#!/usr/bin/env python3
"""
minimap2-Compatible SAM Output

Generates SAM files with all required tags for compatibility with downstream tools
"""

import re
from typing import Dict, List, Optional, Tuple


def parse_cigar_for_nm(cigar: str, query: str, ref: str) -> int:
    """
    Calculate NM tag (edit distance) from CIGAR and sequences
    
    NM = number of mismatches + insertions + deletions
    
    Args:
        cigar: CIGAR string (e.g., "100M2I50M")
        query: Query sequence
        ref: Reference sequence
        
    Returns:
        Edit distance
    """
    if not cigar or cigar == "*":
        return 0
    
    # Parse CIGAR
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    nm = 0
    for length, op in ops:
        length = int(length)
        
        if op in ['I', 'D']:
            # Insertions and deletions count as mismatches
            nm += length
        elif op == 'X':
            # Mismatches
            nm += length
        # M, =, S, H, N, P don't add to NM
    
    return nm


def calculate_mapq(best_score: int, second_score: int, num_seeds: int, 
                   chain_quality: float = 1.0) -> int:
    """
    Calculate MAPQ score (0-60)
    
    Based on minimap2's approach:
    - Unique mappings: MAPQ = 60
    - Ambiguous mappings: MAPQ based on score ratio
    
    Args:
        best_score: Alignment score of best match
        second_score: Alignment score of second-best match (0 if none)
        num_seeds: Number of seeds in chain
        chain_quality: Chain quality metric (0-1)
        
    Returns:
        MAPQ score (0-60)
    """
    if best_score <= 0:
        return 0
    
    # No second-best = unique mapping
    if second_score <= 0:
        return 60
    
    # Calculate score ratio
    score_ratio = best_score / max(second_score, 1)
    
    # Convert to MAPQ
    if score_ratio > 10:
        mapq = 60  # Very confident
    elif score_ratio > 5:
        mapq = 50  # Confident
    elif score_ratio > 3:
        mapq = 40  # Good
    elif score_ratio > 2:
        mapq = 30  # Reasonable
    elif score_ratio > 1.5:
        mapq = 20  # Ambiguous
    elif score_ratio > 1.2:
        mapq = 10  # Very ambiguous
    else:
        mapq = 0   # Multi-mapping
    
    # Adjust for chain quality and seed count
    if num_seeds < 3:
        mapq = max(0, mapq - 10)  # Reduce confidence for few seeds
    
    mapq = int(mapq * chain_quality)  # Scale by chain quality
    
    return min(60, max(0, mapq))


def calculate_alignment_identity(cigar: str) -> float:
    """
    Calculate sequence identity from CIGAR
    
    Identity = matches / (matches + mismatches + insertions + deletions)
    
    Args:
        cigar: CIGAR string
        
    Returns:
        Identity fraction (0-1)
    """
    if not cigar or cigar == "*":
        return 0.0
    
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    matches = 0
    total = 0
    
    for length, op in ops:
        length = int(length)
        
        if op in ['=', 'M']:
            matches += length
            total += length
        elif op in ['X', 'I', 'D']:
            total += length
    
    return matches / total if total > 0 else 0.0


def format_sam_tags(alignment: Dict, second_best_score: int = 0) -> str:
    """
    Generate minimap2-compatible SAM tags
    
    Required tags:
    - NM:i:X - Edit distance
    - AS:i:X - Alignment score
    - ms:i:X - DP alignment score (same as AS for minimap2)
    - nn:i:X - Number of ambiguous bases
    - tp:A:X - Type (P=primary, S=secondary)
    - cm:i:X - Number of minimizer/seeds
    - s1:i:X - Chaining score
    - s2:i:X - Chaining score of second-best
    - de:f:X - Gap-compressed sequence divergence
    - rl:i:X - Length of query regions with repetitive seeds
    
    Args:
        alignment: Alignment dictionary with keys:
            - cigar: CIGAR string
            - alignment_score: Alignment score
            - num_seeds: Number of seeds
            - seed_score: Chain score
            - is_primary: True if primary alignment
        second_best_score: Score of second-best alignment
        
    Returns:
        Tab-separated SAM tags string
    """
    cigar = alignment.get('cigar', '*')
    alignment_score = alignment.get('alignment_score', 0)
    num_seeds = alignment.get('num_seeds', 0)
    seed_score = alignment.get('seed_score', 0)
    is_primary = alignment.get('is_primary', True)
    
    # NM: Edit distance
    nm = alignment.get('nm', 0)
    if nm == 0 and cigar != '*':
        # Calculate from CIGAR if not provided
        nm = parse_cigar_for_nm(cigar, "", "")
    
    # AS/ms: Alignment score
    as_score = int(alignment_score)
    ms_score = as_score  # Same for minimap2
    
    # nn: Ambiguous bases (N's in sequence)
    nn = 0  # Default to 0
    
    # tp: Type (P=primary, S=secondary)
    tp = 'P' if is_primary else 'S'
    
    # cm: Number of seeds/minimizers
    cm = num_seeds
    
    # s1: Chaining score (seed score)
    s1 = int(seed_score)
    
    # s2: Second-best chaining score
    s2 = int(second_best_score)
    
    # de: Divergence
    identity = calculate_alignment_identity(cigar)
    de = 1.0 - identity
    
    # rl: Repetitive seed length (default 0)
    rl = 0
    
    # Format tags
    tags = [
        f"NM:i:{nm}",
        f"ms:i:{ms_score}",
        f"AS:i:{as_score}",
        f"nn:i:{nn}",
        f"tp:A:{tp}",
        f"cm:i:{cm}",
        f"s1:i:{s1}",
        f"s2:i:{s2}",
        f"de:f:{de:.4f}",
        f"rl:i:{rl}"
    ]
    
    return "\t".join(tags)


def format_sam_line(read_id: str, read_seq: str, alignment: Dict, 
                    chr_lengths: Dict[str, int], is_primary: bool = True,
                    second_best_score: int = 0) -> str:
    """
    Format a complete SAM line with all fields and tags
    
    SAM format:
    QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS
    
    Args:
        read_id: Read identifier
        read_seq: Read sequence
        alignment: Alignment dictionary
        chr_lengths: Dict of chromosome lengths
        is_primary: True for primary alignment
        second_best_score: Score of second-best for MAPQ calculation
        
    Returns:
        Complete SAM line
    """
    # Extract alignment info
    chr_name = alignment.get('chr', '*')
    pos = alignment.get('start', 0) + 1  # SAM is 1-based
    cigar = alignment.get('cigar', '*')
    alignment_score = alignment.get('alignment_score', 0)
    num_seeds = alignment.get('num_seeds', 0)
    seed_score = alignment.get('seed_score', 0)
    
    # FLAG
    flag = 0
    if not is_primary:
        flag |= 0x100  # Secondary alignment (256)
    
    # MAPQ
    mapq = calculate_mapq(
        alignment_score,
        second_best_score,
        num_seeds,
        chain_quality=1.0
    )
    
    # RNEXT, PNEXT, TLEN (not used for single-end)
    rnext = '*'
    pnext = 0
    tlen = 0
    
    # QUAL (dummy quality scores)
    qual = 'I' * len(read_seq)  # All high quality
    
    # Prepare alignment dict for tags
    align_with_meta = {
        **alignment,
        'is_primary': is_primary,
        'num_seeds': num_seeds,
        'seed_score': seed_score
    }
    
    # Tags
    tags = format_sam_tags(align_with_meta, second_best_score)
    
    # Assemble SAM line
    sam_line = "\t".join([
        read_id,
        str(flag),
        chr_name,
        str(pos),
        str(mapq),
        cigar,
        rnext,
        str(pnext),
        str(tlen),
        read_seq,
        qual,
        tags
    ])
    
    return sam_line


def write_sam_header(chr_lengths: Dict[str, int]) -> List[str]:
    """
    Generate SAM header lines
    
    Args:
        chr_lengths: Dict of chromosome names → lengths
        
    Returns:
        List of header lines
    """
    header = []
    
    # @HD line
    header.append("@HD\tVN:1.6\tSO:unsorted")
    
    # @SQ lines (sorted for consistency)
    for chr_name in sorted(chr_lengths.keys()):
        length = chr_lengths[chr_name]
        header.append(f"@SQ\tSN:{chr_name}\tLN:{length}")
    
    # @PG line
    header.append("@PG\tID:genocache\tPN:genocache\tVN:4.1\tCL:genocache_align.py")
    
    return header


def write_sam_with_secondaries(output_path: str, alignments: List[Dict],
                               chr_lengths: Dict[str, int]):
    """
    Write SAM file with primary and secondary alignments
    
    Args:
        output_path: Output SAM file path
        alignments: List of alignment dicts, each with:
            - read_id: Read identifier
            - read_seq: Read sequence
            - results: List of alignment results (sorted by score)
        chr_lengths: Dict of chromosome lengths
    """
    with open(output_path, 'w') as f:
        # Write header
        for line in write_sam_header(chr_lengths):
            f.write(line + "\n")
        
        # Write alignments
        for aln in alignments:
            read_id = aln['read_id']
            read_seq = aln['read_seq']
            results = aln['results']  # Sorted by score (best first)
            
            if not results:
                # Unmapped
                f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t{'I'*len(read_seq)}\n")
                continue
            
            # Primary alignment (best)
            primary = results[0]
            second_score = results[1]['alignment_score'] if len(results) > 1 else 0
            
            sam_line = format_sam_line(
                read_id, read_seq, primary, chr_lengths,
                is_primary=True, second_best_score=second_score
            )
            f.write(sam_line + "\n")
            
            # Secondary alignments (if any)
            for secondary in results[1:]:
                # Only output if score is reasonable
                if secondary['alignment_score'] > primary['alignment_score'] * 0.5:
                    sam_line = format_sam_line(
                        read_id, read_seq, secondary, chr_lengths,
                        is_primary=False, second_best_score=0
                    )
                    f.write(sam_line + "\n")


# Test function
if __name__ == '__main__':
    print("SAM Output Module - minimap2 Compatible")
    print("="*60)
    print()
    
    # Test tag generation
    test_alignment = {
        'cigar': '100M2I50M',
        'alignment_score': 950,
        'num_seeds': 5,
        'seed_score': 25,
        'is_primary': True
    }
    
    tags = format_sam_tags(test_alignment, second_best_score=800)
    print("Test tags:")
    print(tags)
    print()
    
    # Test MAPQ calculation
    mapq1 = calculate_mapq(1000, 100, 5)  # Clear winner
    mapq2 = calculate_mapq(1000, 900, 5)  # Ambiguous
    mapq3 = calculate_mapq(1000, 0, 5)    # Unique
    
    print(f"MAPQ (score 1000 vs 100, 5 seeds): {mapq1}")
    print(f"MAPQ (score 1000 vs 900, 5 seeds): {mapq2}")
    print(f"MAPQ (score 1000 vs 0, 5 seeds):   {mapq3}")
    print()
    
    print("✅ SAM output module ready!")
