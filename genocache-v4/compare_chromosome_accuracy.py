#!/usr/bin/env python3
"""
Compare Chromosome Accuracy - Simple Analysis

This script compares GenoCache vs minimap2 chromosome-level accuracy
using existing SAM files (no model loading needed).

Expected: Shows our 37% accuracy bug
After EXTEND fix: Should show 95%+ accuracy
"""

import sys
from pathlib import Path

def parse_sam_alignments(sam_file):
    """Parse SAM file and extract alignments"""
    alignments = {}
    
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            read_id = fields[0]
            flag = int(fields[1])
            ref_chr = fields[2]
            pos = int(fields[3])
            
            # Skip unmapped
            if flag & 4:
                continue
            
            alignments[read_id] = {
                'chr': ref_chr,
                'pos': pos,
                'flag': flag
            }
    
    return alignments

def compare_aligners(genocache_sam, minimap2_sam):
    """Compare chromosome-level accuracy"""
    print("=" * 80)
    print("CHROMOSOME ACCURACY COMPARISON")
    print("=" * 80)
    print()
    
    # Parse both SAM files
    print(f"Parsing GenoCache SAM: {genocache_sam}")
    genocache = parse_sam_alignments(genocache_sam)
    print(f"  → {len(genocache)} mapped reads")
    
    print(f"Parsing minimap2 SAM: {minimap2_sam}")
    minimap2 = parse_sam_alignments(minimap2_sam)
    print(f"  → {len(minimap2)} mapped reads")
    print()
    
    # Find common reads
    common_reads = set(genocache.keys()) & set(minimap2.keys())
    print(f"Common mapped reads: {len(common_reads)}")
    print()
    
    if not common_reads:
        print("❌ No common reads to compare!")
        return
    
    # Compare chromosomes
    chr_match = 0
    chr_mismatch = 0
    
    print("DETAILED COMPARISON:")
    print("-" * 80)
    
    for read_id in sorted(common_reads):
        gc_chr = genocache[read_id]['chr']
        mm_chr = minimap2[read_id]['chr']
        gc_pos = genocache[read_id]['pos']
        mm_pos = minimap2[read_id]['pos']
        
        if gc_chr == mm_chr:
            chr_match += 1
            pos_diff = abs(gc_pos - mm_pos)
            status = "✅ MATCH"
            print(f"{read_id}:")
            print(f"  GenoCache:  {gc_chr}:{gc_pos:,}")
            print(f"  minimap2:   {mm_chr}:{mm_pos:,}")
            print(f"  {status} (position diff: {pos_diff:,} bp)")
        else:
            chr_mismatch += 1
            status = "❌ CHR MISMATCH"
            print(f"{read_id}:")
            print(f"  GenoCache:  {gc_chr}:{gc_pos:,}")
            print(f"  minimap2:   {mm_chr}:{mm_pos:,}")
            print(f"  {status}")
        print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    total = len(common_reads)
    chr_accuracy = chr_match / total * 100
    
    print(f"Total compared: {total}")
    print(f"Chromosome matches: {chr_match} ({chr_accuracy:.1f}%)")
    print(f"Chromosome mismatches: {chr_mismatch} ({chr_mismatch/total*100:.1f}%)")
    print()
    
    if chr_accuracy < 50:
        print("⚠️  CRITICAL: Chromosome accuracy < 50%!")
        print("    This indicates a serious bug in chromosome selection.")
        print("    → Need to implement EXTEND phase (align to all candidates)")
    elif chr_accuracy < 90:
        print("⚠️  WARNING: Chromosome accuracy < 90%")
        print("    Acceptable for research but needs improvement for production.")
    else:
        print("✅ GOOD: Chromosome accuracy > 90%")
        print("    Quality is acceptable for production use.")
    
    print()
    
    # Position accuracy (for correctly mapped reads)
    if chr_match > 0:
        print("POSITION ACCURACY (for correct chromosomes):")
        print("-" * 80)
        
        pos_diffs = []
        for read_id in common_reads:
            if genocache[read_id]['chr'] == minimap2[read_id]['chr']:
                diff = abs(genocache[read_id]['pos'] - minimap2[read_id]['pos'])
                pos_diffs.append(diff)
        
        pos_diffs.sort()
        median_diff = pos_diffs[len(pos_diffs)//2]
        
        print(f"Median position difference: {median_diff:,} bp")
        print(f"Min: {min(pos_diffs):,} bp")
        print(f"Max: {max(pos_diffs):,} bp")
        print()
        
        within_1kb = sum(1 for d in pos_diffs if d <= 1000)
        print(f"Within 1kb: {within_1kb}/{len(pos_diffs)} ({within_1kb/len(pos_diffs)*100:.1f}%)")
    
    return {
        'chr_accuracy': chr_accuracy,
        'chr_match': chr_match,
        'chr_mismatch': chr_mismatch,
        'total': total
    }

if __name__ == '__main__':
    # Check for existing SAM files
    genocache_sam = Path("test_complete_10reads.sam")
    minimap2_sam = Path("minimap2_same_10reads.sam")
    
    if not genocache_sam.exists():
        print(f"❌ GenoCache SAM not found: {genocache_sam}")
        print("   Need to run the complete pipeline first!")
        sys.exit(1)
    
    if not minimap2_sam.exists():
        print(f"❌ minimap2 SAM not found: {minimap2_sam}")
        print("   Need to run minimap2 on same reads!")
        sys.exit(1)
    
    compare_aligners(genocache_sam, minimap2_sam)
