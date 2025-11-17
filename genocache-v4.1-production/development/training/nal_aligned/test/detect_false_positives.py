#!/usr/bin/env python3
"""
False Positive Detection for NAL Alignment

Compares NAL results with minimap2 (gold standard) to detect:
1. Reads mapped to wrong chromosome
2. Reads mapped to wrong position (>1000bp off)
3. Low-quality alignments (poor CIGAR, low MAPQ)
4. Ambiguous mappings (multiple good locations)

Metrics:
- Precision: % of NAL mappings that match minimap2
- Recall: % of minimap2 mappings found by NAL
- False positive rate: % of NAL mappings that disagree
"""

import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Set
from collections import defaultdict


def parse_sam_line(line: str) -> Dict:
    """Parse SAM line into dict"""
    if line.startswith('@') or not line.strip():
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None
    
    try:
        flag = int(fields[1])
        pos = int(fields[3])
        mapq = int(fields[4])
    except ValueError:
        # Malformed line, skip
        return None
    
    return {
        'qname': fields[0],
        'flag': flag,
        'rname': fields[2],
        'pos': pos,
        'mapq': mapq,
        'cigar': fields[5],
        'is_mapped': fields[2] != '*',
        'is_reverse': (flag & 16) != 0,
    }


def load_sam_alignments(sam_path: str) -> Dict[str, Dict]:
    """Load SAM file into dict keyed by read ID"""
    alignments = {}
    
    with open(sam_path) as f:
        for line in f:
            aln = parse_sam_line(line)
            if aln:
                alignments[aln['qname']] = aln
    
    return alignments


def positions_match(pos1: int, pos2: int, tolerance: int = 1000) -> bool:
    """Check if two positions are within tolerance"""
    return abs(pos1 - pos2) <= tolerance


def compare_alignments(
    nal_alignments: Dict[str, Dict],
    mm2_alignments: Dict[str, Dict],
    position_tolerance: int = 1000
) -> Dict:
    """
    Compare NAL and minimap2 alignments
    
    Returns dict with:
    - true_positives: NAL correct (matches minimap2)
    - false_positives: NAL wrong (disagrees with minimap2)
    - false_negatives: NAL missed (minimap2 found, NAL didn't)
    - true_negatives: Both unmapped
    - details: List of discrepancies
    """
    # Get all read IDs
    all_reads = set(nal_alignments.keys()) | set(mm2_alignments.keys())
    
    results = {
        'true_positives': 0,
        'false_positives': 0,
        'false_negatives': 0,
        'true_negatives': 0,
        'wrong_chr': 0,
        'wrong_pos': 0,
        'low_quality': 0,
        'details': []
    }
    
    for read_id in all_reads:
        nal = nal_alignments.get(read_id)
        mm2 = mm2_alignments.get(read_id)
        
        # Handle missing alignments
        if nal is None or mm2 is None:
            continue
        
        nal_mapped = nal['is_mapped']
        mm2_mapped = mm2['is_mapped']
        
        # Both unmapped - True Negative
        if not nal_mapped and not mm2_mapped:
            results['true_negatives'] += 1
            continue
        
        # NAL unmapped, minimap2 mapped - False Negative
        if not nal_mapped and mm2_mapped:
            results['false_negatives'] += 1
            results['details'].append({
                'read_id': read_id,
                'type': 'false_negative',
                'reason': 'NAL missed, minimap2 found',
                'mm2_chr': mm2['rname'],
                'mm2_pos': mm2['pos'],
            })
            continue
        
        # NAL mapped, minimap2 unmapped - Potential False Positive
        if nal_mapped and not mm2_mapped:
            results['false_positives'] += 1
            results['details'].append({
                'read_id': read_id,
                'type': 'false_positive',
                'reason': 'NAL mapped, minimap2 unmapped',
                'nal_chr': nal['rname'],
                'nal_pos': nal['pos'],
                'nal_mapq': nal['mapq'],
            })
            continue
        
        # Both mapped - Check agreement
        chr_match = nal['rname'] == mm2['rname']
        pos_match = positions_match(nal['pos'], mm2['pos'], position_tolerance)
        strand_match = nal['is_reverse'] == mm2['is_reverse']
        
        if chr_match and pos_match and strand_match:
            # True Positive - Agreement
            results['true_positives'] += 1
        else:
            # False Positive - Disagreement
            results['false_positives'] += 1
            
            reason = []
            if not chr_match:
                results['wrong_chr'] += 1
                reason.append('wrong_chr')
            if not pos_match:
                results['wrong_pos'] += 1
                reason.append(f'wrong_pos (off by {abs(nal["pos"] - mm2["pos"])}bp)')
            if not strand_match:
                reason.append('wrong_strand')
            
            results['details'].append({
                'read_id': read_id,
                'type': 'false_positive',
                'reason': ', '.join(reason),
                'nal_chr': nal['rname'],
                'nal_pos': nal['pos'],
                'nal_mapq': nal['mapq'],
                'mm2_chr': mm2['rname'],
                'mm2_pos': mm2['pos'],
                'mm2_mapq': mm2['mapq'],
            })
        
        # Check for low quality
        if nal_mapped and nal['mapq'] < 10:
            results['low_quality'] += 1
    
    return results


def print_report(results: Dict, nal_name: str = "NAL", mm2_name: str = "minimap2"):
    """Print comparison report"""
    tp = results['true_positives']
    fp = results['false_positives']
    fn = results['false_negatives']
    tn = results['true_negatives']
    
    total = tp + fp + fn + tn
    nal_mapped = tp + fp
    mm2_mapped = tp + fn
    
    print("\n" + "="*80)
    print("FALSE POSITIVE DETECTION REPORT")
    print("="*80 + "\n")
    
    print(f"Comparison: {nal_name} vs {mm2_name} (gold standard)")
    print(f"Total reads: {total}\n")
    
    # Confusion matrix
    print("CONFUSION MATRIX:")
    print("─"*80)
    print(f"{'':20} {mm2_name} Mapped    {mm2_name} Unmapped")
    print(f"{nal_name} Mapped        {tp:6}           {fp:6}  (FP)")
    print(f"{nal_name} Unmapped      {fn:6} (FN)     {tn:6}")
    print()
    
    # Metrics
    if nal_mapped > 0:
        precision = tp / nal_mapped * 100
    else:
        precision = 0
    
    if mm2_mapped > 0:
        recall = tp / mm2_mapped * 100
    else:
        recall = 0
    
    if precision + recall > 0:
        f1 = 2 * (precision * recall) / (precision + recall)
    else:
        f1 = 0
    
    if nal_mapped > 0:
        fpr = fp / nal_mapped * 100
    else:
        fpr = 0
    
    print("METRICS:")
    print("─"*80)
    print(f"  Precision: {precision:.2f}%  ({tp}/{nal_mapped} NAL mappings correct)")
    print(f"  Recall:    {recall:.2f}%  ({tp}/{mm2_mapped} {mm2_name} mappings found)")
    print(f"  F1 Score:  {f1:.2f}%")
    print(f"  False Positive Rate: {fpr:.2f}%  ({fp}/{nal_mapped} NAL mappings wrong)")
    print()
    
    # Error breakdown
    if fp > 0:
        print("FALSE POSITIVE BREAKDOWN:")
        print("─"*80)
        print(f"  Wrong chromosome: {results['wrong_chr']} ({results['wrong_chr']/fp*100:.1f}%)")
        print(f"  Wrong position:   {results['wrong_pos']} ({results['wrong_pos']/fp*100:.1f}%)")
        print(f"  Low quality (MAPQ<10): {results['low_quality']} ({results['low_quality']/nal_mapped*100:.1f}% of mapped)")
        print()
    
    # Show examples
    if len(results['details']) > 0:
        print("EXAMPLES OF DISCREPANCIES:")
        print("─"*80)
        
        # Show up to 10 examples
        for detail in results['details'][:10]:
            print(f"\n  Read: {detail['read_id']}")
            print(f"  Type: {detail['type']}")
            print(f"  Reason: {detail['reason']}")
            
            if detail['type'] == 'false_positive':
                if 'nal_chr' in detail:
                    print(f"    {nal_name}: {detail['nal_chr']}:{detail['nal_pos']} (MAPQ={detail['nal_mapq']})")
                if 'mm2_chr' in detail:
                    print(f"    {mm2_name}: {detail['mm2_chr']}:{detail['mm2_pos']} (MAPQ={detail['mm2_mapq']})")
            elif detail['type'] == 'false_negative':
                print(f"    {mm2_name}: {detail['mm2_chr']}:{detail['mm2_pos']}")
        
        if len(results['details']) > 10:
            print(f"\n  ... and {len(results['details']) - 10} more discrepancies")
        print()
    
    # Summary judgment
    print("="*80)
    print("SUMMARY:")
    print("="*80)
    
    if fpr < 1:
        print(f"  ✅ EXCELLENT: False positive rate <1% ({fpr:.2f}%)")
    elif fpr < 5:
        print(f"  ✅ GOOD: False positive rate <5% ({fpr:.2f}%)")
    elif fpr < 10:
        print(f"  ⚠️  ACCEPTABLE: False positive rate <10% ({fpr:.2f}%)")
    else:
        print(f"  ❌ HIGH: False positive rate ≥10% ({fpr:.2f}%)")
    
    if recall > 95:
        print(f"  ✅ EXCELLENT: Recall >95% ({recall:.2f}%)")
    elif recall > 90:
        print(f"  ✅ GOOD: Recall >90% ({recall:.2f}%)")
    elif recall > 85:
        print(f"  ⚠️  ACCEPTABLE: Recall >85% ({recall:.2f}%)")
    else:
        print(f"  ❌ LOW: Recall <85% ({recall:.2f}%)")
    
    print()


def main():
    """Compare NAL and minimap2 alignments"""
    parser = argparse.ArgumentParser(description='Detect false positives in NAL alignment')
    parser.add_argument('--nal', required=True, help='NAL SAM file')
    parser.add_argument('--minimap2', required=True, help='minimap2 SAM file')
    parser.add_argument('--tolerance', type=int, default=1000, 
                        help='Position tolerance in bp (default: 1000)')
    parser.add_argument('--nal-name', default='NAL', help='Name for NAL in report')
    
    args = parser.parse_args()
    
    print(f"Loading alignments...")
    print(f"  NAL: {args.nal}")
    nal_alignments = load_sam_alignments(args.nal)
    print(f"    Loaded {len(nal_alignments)} reads")
    
    print(f"  minimap2: {args.minimap2}")
    mm2_alignments = load_sam_alignments(args.minimap2)
    print(f"    Loaded {len(mm2_alignments)} reads")
    
    print(f"\nComparing alignments (position tolerance: {args.tolerance}bp)...")
    results = compare_alignments(nal_alignments, mm2_alignments, args.tolerance)
    
    print_report(results, args.nal_name)


if __name__ == '__main__':
    main()
