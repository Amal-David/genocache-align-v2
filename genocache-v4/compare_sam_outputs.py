#!/usr/bin/env python3
"""
Head-to-head comparison of GenoCache vs minimap2 SAM outputs
"""

import re
from collections import defaultdict

def parse_sam_line(line):
    """Parse SAM line to extract key fields"""
    if line.startswith('@'):
        return None  # Header
    
    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None
    
    return {
        'qname': fields[0],
        'flag': int(fields[1]),
        'rname': fields[2],
        'pos': int(fields[3]),
        'mapq': int(fields[4]),
        'cigar': fields[5],
        'seq': fields[9],
        'is_mapped': fields[2] != '*'
    }

def parse_cigar(cigar):
    """Parse CIGAR to get operation counts"""
    if cigar == '*':
        return {}
    
    ops = defaultdict(int)
    pattern = r'(\d+)([MIDNSHPX=])'
    for match in re.finditer(pattern, cigar):
        count = int(match.group(1))
        op = match.group(2)
        ops[op] += count
    
    return dict(ops)

def compare_alignments(geocache_sam, minimap2_sam):
    """Compare two SAM files line by line"""
    
    print("=" * 80)
    print("GenoCache vs minimap2 - Line-by-Line Comparison")
    print("=" * 80)
    print()
    
    # Parse GenoCache
    geocache = {}
    with open(geocache_sam) as f:
        for line in f:
            rec = parse_sam_line(line)
            if rec:
                geocache[rec['qname']] = rec
    
    # Parse minimap2
    minimap2 = {}
    with open(minimap2_sam) as f:
        for line in f:
            rec = parse_sam_line(line)
            if rec:
                minimap2[rec['qname']] = rec
    
    print(f"GenoCache: {len(geocache)} alignments")
    print(f"minimap2:  {len(minimap2)} alignments")
    print()
    
    # Compare each read
    all_reads = set(geocache.keys()) | set(minimap2.keys())
    
    stats = {
        'both_mapped': 0,
        'only_geocache': 0,
        'only_minimap2': 0,
        'both_unmapped': 0,
        'same_chr': 0,
        'diff_chr': 0,
        'pos_diff': []
    }
    
    print("=" * 80)
    print("DETAILED COMPARISON")
    print("=" * 80)
    print()
    
    for read_id in sorted(all_reads)[:10]:  # First 10 reads
        print(f"Read: {read_id}")
        print("-" * 80)
        
        gc = geocache.get(read_id)
        mm = minimap2.get(read_id)
        
        if not gc and not mm:
            print("  Both: UNMAPPED")
            stats['both_unmapped'] += 1
        elif not gc:
            print(f"  GenoCache: UNMAPPED")
            print(f"  minimap2:  {mm['rname']}:{mm['pos']} (CIGAR: {mm['cigar'][:50]}...)")
            stats['only_minimap2'] += 1
        elif not mm:
            print(f"  GenoCache: {gc['rname']}:{gc['pos']} (CIGAR: {gc['cigar'][:50]}...)")
            print(f"  minimap2:  UNMAPPED")
            stats['only_geocache'] += 1
        else:
            # Both mapped - compare details
            if not gc['is_mapped']:
                print(f"  GenoCache: UNMAPPED")
            else:
                print(f"  GenoCache: {gc['rname']}:{gc['pos']:,} (MAPQ: {gc['mapq']})")
                print(f"    CIGAR: {gc['cigar'][:70]}{'...' if len(gc['cigar']) > 70 else ''}")
                gc_ops = parse_cigar(gc['cigar'])
                print(f"    Ops: M={gc_ops.get('M',0)}, I={gc_ops.get('I',0)}, D={gc_ops.get('D',0)}")
            
            if not mm['is_mapped']:
                print(f"  minimap2:  UNMAPPED")
            else:
                print(f"  minimap2:  {mm['rname']}:{mm['pos']:,} (MAPQ: {mm['mapq']})")
                print(f"    CIGAR: {mm['cigar'][:70]}{'...' if len(mm['cigar']) > 70 else ''}")
                mm_ops = parse_cigar(mm['cigar'])
                print(f"    Ops: M={mm_ops.get('M',0)}, I={mm_ops.get('I',0)}, D={mm_ops.get('D',0)}")
            
            if gc['is_mapped'] and mm['is_mapped']:
                stats['both_mapped'] += 1
                
                # Compare chromosome
                if gc['rname'] == mm['rname']:
                    stats['same_chr'] += 1
                    print(f"  ✅ Same chromosome: {gc['rname']}")
                    
                    # Compare position
                    pos_diff = abs(gc['pos'] - mm['pos'])
                    stats['pos_diff'].append(pos_diff)
                    print(f"  Position diff: {pos_diff:,} bp")
                    
                    if pos_diff < 1000:
                        print(f"  ✅ Close alignment (within 1kb)")
                    else:
                        print(f"  ⚠️  Distant alignment (>{pos_diff/1000:.1f}kb apart)")
                else:
                    stats['diff_chr'] += 1
                    print(f"  ❌ Different chromosomes!")
                    print(f"     GenoCache: {gc['rname']}")
                    print(f"     minimap2:  {mm['rname']}")
        
        print()
    
    # Summary statistics
    print("=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print()
    
    total = len(all_reads)
    print(f"Total reads: {total}")
    print()
    
    print("Mapping status:")
    print(f"  Both mapped:     {stats['both_mapped']:3d} ({100*stats['both_mapped']/total:.1f}%)")
    print(f"  Only GenoCache:  {stats['only_geocache']:3d} ({100*stats['only_geocache']/total:.1f}%)")
    print(f"  Only minimap2:   {stats['only_minimap2']:3d} ({100*stats['only_minimap2']/total:.1f}%)")
    print(f"  Both unmapped:   {stats['both_unmapped']:3d} ({100*stats['both_unmapped']/total:.1f}%)")
    print()
    
    if stats['both_mapped'] > 0:
        print("For reads mapped by both:")
        print(f"  Same chromosome: {stats['same_chr']}/{stats['both_mapped']} ({100*stats['same_chr']/stats['both_mapped']:.1f}%)")
        print(f"  Diff chromosome: {stats['diff_chr']}/{stats['both_mapped']} ({100*stats['diff_chr']/stats['both_mapped']:.1f}%)")
        
        if stats['pos_diff']:
            print()
            print("Position differences:")
            print(f"  Min:    {min(stats['pos_diff']):,} bp")
            print(f"  Max:    {max(stats['pos_diff']):,} bp")
            print(f"  Median: {sorted(stats['pos_diff'])[len(stats['pos_diff'])//2]:,} bp")
            print(f"  Mean:   {sum(stats['pos_diff'])/len(stats['pos_diff']):,.0f} bp")
    
    print()
    print("=" * 80)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 2:
        compare_alignments(sys.argv[1], sys.argv[2])
    else:
        compare_alignments('test_complete_10reads.sam', 'minimap2_full.sam')
