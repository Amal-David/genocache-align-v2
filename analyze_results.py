#!/usr/bin/env python3
"""
Analyze alignment results and generate statistics
"""
import sys
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

def load_alignments(tsv_file):
    """Load alignment results"""
    print(f"Loading alignments from {tsv_file}...")
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"  ✓ Loaded {len(df):,} alignments")
    print(f"  ✓ Unique reads: {df['read_id'].nunique():,}")
    return df

def analyze_top1(df):
    """Analyze top-1 alignments"""
    print("\n" + "="*70)
    print("Top-1 Alignment Analysis")
    print("="*70)
    
    top1 = df[df['rank'] == 1].copy()
    
    print(f"\nTotal reads with top-1 alignment: {len(top1):,}")
    
    # Score statistics
    print(f"\nScore statistics:")
    print(f"  Mean:   {top1['score'].mean():.4f}")
    print(f"  Median: {top1['score'].median():.4f}")
    print(f"  Min:    {top1['score'].min():.4f}")
    print(f"  Max:    {top1['score'].max():.4f}")
    print(f"  Std:    {top1['score'].std():.4f}")
    
    # Score distribution
    high_conf = len(top1[top1['score'] > 0.8])
    med_conf = len(top1[(top1['score'] > 0.6) & (top1['score'] <= 0.8)])
    low_conf = len(top1[top1['score'] <= 0.6])
    
    print(f"\nConfidence distribution:")
    print(f"  High (>0.8):  {high_conf:,} ({high_conf/len(top1)*100:.1f}%)")
    print(f"  Medium (0.6-0.8): {med_conf:,} ({med_conf/len(top1)*100:.1f}%)")
    print(f"  Low (<0.6):   {low_conf:,} ({low_conf/len(top1)*100:.1f}%)")
    
    return top1

def analyze_chromosomes(top1):
    """Analyze chromosome distribution"""
    print(f"\n" + "="*70)
    print("Chromosome Distribution")
    print("="*70)
    
    chr_counts = top1['chrom'].value_counts().sort_index()
    
    print(f"\nReads per chromosome:")
    print(f"{'Chrom':>6} {'Count':>10} {'Percent':>8}")
    print("-" * 26)
    
    total = len(top1)
    for chrom, count in chr_counts.items():
        pct = count / total * 100
        print(f"{chrom:>6} {count:>10,} {pct:>7.2f}%")
    
    # Expected vs actual (rough estimate based on chromosome sizes)
    # GRCh38 chromosome sizes (approximate, Mb)
    chr_sizes = {
        0: 249, 1: 242, 2: 198, 3: 191, 4: 181, 5: 171,
        6: 171, 7: 159, 8: 146, 9: 141, 10: 136, 11: 135,
        12: 133, 13: 115, 14: 107, 15: 102, 16: 90, 17: 83,
        18: 80, 19: 59, 20: 64, 21: 47, 22: 51, 23: 156, 24: 57
    }
    
    print(f"\nExpected vs Actual (top 10 chromosomes):")
    print(f"{'Chrom':>6} {'Expected%':>10} {'Actual%':>10} {'Ratio':>8}")
    print("-" * 36)
    
    total_size = sum(chr_sizes.values())
    for chrom in sorted(chr_counts.index)[:10]:
        if chrom in chr_sizes:
            expected_pct = chr_sizes[chrom] / total_size * 100
            actual_pct = chr_counts[chrom] / total * 100
            ratio = actual_pct / expected_pct if expected_pct > 0 else 0
            print(f"{chrom:>6} {expected_pct:>9.2f}% {actual_pct:>9.2f}% {ratio:>7.2f}x")

def analyze_positions(top1):
    """Analyze position distribution"""
    print(f"\n" + "="*70)
    print("Position Distribution")
    print("="*70)
    
    print(f"\nPosition statistics:")
    print(f"  Min:    {top1['position'].min():,}")
    print(f"  Max:    {top1['position'].max():,}")
    print(f"  Mean:   {top1['position'].mean():,.0f}")
    print(f"  Median: {top1['position'].median():,.0f}")
    
    # Position coverage per chromosome
    print(f"\nPosition range per chromosome (top 10):")
    print(f"{'Chrom':>6} {'Min Pos':>12} {'Max Pos':>12} {'Range (Mb)':>12}")
    print("-" * 44)
    
    for chrom in sorted(top1['chrom'].unique())[:10]:
        chr_data = top1[top1['chrom'] == chrom]
        min_pos = chr_data['position'].min()
        max_pos = chr_data['position'].max()
        range_mb = (max_pos - min_pos) / 1e6
        print(f"{chrom:>6} {min_pos:>12,} {max_pos:>12,} {range_mb:>11.1f}")

def analyze_multi_mapping(df):
    """Analyze reads with multiple good alignments"""
    print(f"\n" + "="*70)
    print("Multi-mapping Analysis")
    print("="*70)
    
    # Group by read
    read_groups = df.groupby('read_id')
    
    # High-confidence alignments per read (score > 0.7)
    high_conf_counts = []
    for read_id, group in read_groups:
        high_conf = len(group[group['score'] > 0.7])
        high_conf_counts.append(high_conf)
    
    high_conf_counts = np.array(high_conf_counts)
    
    print(f"\nHigh-confidence alignments per read (score > 0.7):")
    print(f"  Mean:   {high_conf_counts.mean():.2f}")
    print(f"  Median: {np.median(high_conf_counts):.0f}")
    print(f"  Max:    {high_conf_counts.max()}")
    
    # Distribution
    unique_counts = Counter(high_conf_counts)
    print(f"\nDistribution:")
    for count in sorted(unique_counts.keys())[:10]:
        n_reads = unique_counts[count]
        pct = n_reads / len(high_conf_counts) * 100
        print(f"  {count} alignment(s): {n_reads:,} reads ({pct:.1f}%)")
    
    # Ambiguous reads (multiple high-confidence alignments)
    ambiguous = len(high_conf_counts[high_conf_counts > 1])
    print(f"\nAmbiguous reads (>1 high-conf): {ambiguous:,} ({ambiguous/len(high_conf_counts)*100:.1f}%)")

def generate_summary_report(df, output_file='analysis_summary.txt'):
    """Generate text summary report"""
    print(f"\n" + "="*70)
    print("Generating Summary Report")
    print("="*70)
    
    with open(output_file, 'w') as f:
        f.write("GenoCache Alignment Analysis Report\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total alignments: {len(df):,}\n")
        f.write(f"Unique reads: {df['read_id'].nunique():,}\n")
        f.write(f"Alignments per read: {len(df) / df['read_id'].nunique():.1f}\n\n")
        
        top1 = df[df['rank'] == 1]
        f.write(f"Top-1 Statistics:\n")
        f.write(f"  Mean score: {top1['score'].mean():.4f}\n")
        f.write(f"  Median score: {top1['score'].median():.4f}\n")
        f.write(f"  High confidence (>0.8): {len(top1[top1['score'] > 0.8]):,}\n")
        f.write(f"  Unique chromosomes: {top1['chrom'].nunique()}\n\n")
        
        f.write("Top 10 chromosomes:\n")
        chr_counts = top1['chrom'].value_counts()
        for i, (chrom, count) in enumerate(chr_counts.head(10).items()):
            f.write(f"  {i+1}. Chr {chrom}: {count:,} reads\n")
    
    print(f"  ✓ Saved to {output_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyze GenoCache alignment results')
    parser.add_argument('--input', required=True, help='Alignments TSV file')
    parser.add_argument('--output', default='analysis_summary.txt', help='Summary output file')
    args = parser.parse_args()
    
    # Load data
    df = load_alignments(args.input)
    
    # Run analyses
    top1 = analyze_top1(df)
    analyze_chromosomes(top1)
    analyze_positions(top1)
    analyze_multi_mapping(df)
    
    # Generate report
    generate_summary_report(df, args.output)
    
    print(f"\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)

if __name__ == '__main__':
    main()
