#!/usr/bin/env python3
"""
Analyze Speed vs Accuracy Trade-offs

Finds optimal configurations for different use cases:
- Maximum accuracy (regardless of speed)
- Best speed/accuracy balance
- Fastest while maintaining >90% accuracy
"""

import csv
import sys
from pathlib import Path
from typing import List, Dict
import numpy as np


def load_csv(filepath: str) -> List[Dict]:
    """Load CSV results"""
    results = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for key in row:
                try:
                    if '.' in str(row[key]):
                        row[key] = float(row[key])
                    else:
                        row[key] = int(row[key])
                except:
                    pass
            results.append(row)
    return results


def analyze_speed_accuracy(results: List[Dict]):
    """Analyze speed vs accuracy trade-offs"""
    print("\n" + "="*80)
    print("SPEED vs ACCURACY ANALYSIS")
    print("="*80 + "\n")
    
    # Sort by seeds, then K
    results = sorted(results, key=lambda x: (x['num_seeds'], x['K']))
    
    # Display full table
    print(f"{'Seeds':<8} {'K':<6} {'Mapping':<10} {'Time':<10} {'Enc':<8} {'FAISS':<8} {'Chain':<8} {'Anchors':<10} {'Chain':<8}")
    print(f"{'':8} {'':6} {'%':10} {'(ms)':10} {'(ms)':8} {'(ms)':8} {'(ms)':8} {'/seed':10} {'Score':8}")
    print("-" * 100)
    
    for r in results:
        print(f"{r['num_seeds']:<8} {r['K']:<6} {r['mapping_rate']:<10.1f} "
              f"{r['total_time_ms']:<10.1f} {r['encoding_time_ms']:<8.1f} "
              f"{r['faiss_time_ms']:<8.1f} {r['chaining_time_ms']:<8.1f} "
              f"{r['avg_anchors_per_seed']:<10.1f} {r['avg_best_chain_score']:<8.1f}")
    
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80 + "\n")
    
    # Find key configurations
    best_accuracy = max(results, key=lambda x: x['mapping_rate'])
    fastest = min(results, key=lambda x: x['total_time_ms'])
    
    # Best at each K value
    print("📊 Best mapping at each K:")
    for K in sorted(set(r['K'] for r in results)):
        k_results = [r for r in results if r['K'] == K]
        best_k = max(k_results, key=lambda x: x['mapping_rate'])
        print(f"   K={K}: seeds={best_k['num_seeds']} → {best_k['mapping_rate']:.1f}% "
              f"({best_k['total_time_ms']:.1f}ms)")
    
    print("\n📊 Best mapping at each seed count:")
    for seeds in sorted(set(r['num_seeds'] for r in results)):
        seed_results = [r for r in results if r['num_seeds'] == seeds]
        best_seed = max(seed_results, key=lambda x: x['mapping_rate'])
        print(f"   seeds={seeds}: K={best_seed['K']} → {best_seed['mapping_rate']:.1f}% "
              f"({best_seed['total_time_ms']:.1f}ms)")
    
    # Speed/accuracy trade-offs
    print("\n" + "="*80)
    print("RECOMMENDED CONFIGURATIONS")
    print("="*80 + "\n")
    
    # Find configs with >90% accuracy
    good_accuracy = [r for r in results if r['mapping_rate'] >= 90]
    
    if good_accuracy:
        # Best speed with good accuracy
        fastest_good = min(good_accuracy, key=lambda x: x['total_time_ms'])
        
        print("🚀 FASTEST with ≥90% accuracy:")
        print(f"   Seeds: {fastest_good['num_seeds']}, K: {fastest_good['K']}")
        print(f"   Mapping: {fastest_good['mapping_rate']:.1f}%")
        print(f"   Time: {fastest_good['total_time_ms']:.1f}ms/read")
        print(f"   Breakdown: enc={fastest_good['encoding_time_ms']:.1f}ms, "
              f"faiss={fastest_good['faiss_time_ms']:.1f}ms, "
              f"chain={fastest_good['chaining_time_ms']:.1f}ms")
        print()
    
    # Best accuracy
    print("🎯 MAXIMUM ACCURACY:")
    print(f"   Seeds: {best_accuracy['num_seeds']}, K: {best_accuracy['K']}")
    print(f"   Mapping: {best_accuracy['mapping_rate']:.1f}%")
    print(f"   Time: {best_accuracy['total_time_ms']:.1f}ms/read")
    print(f"   Breakdown: enc={best_accuracy['encoding_time_ms']:.1f}ms, "
          f"faiss={best_accuracy['faiss_time_ms']:.1f}ms, "
          f"chain={best_accuracy['chaining_time_ms']:.1f}ms")
    print()
    
    # Find sweet spot (efficiency metric)
    for r in results:
        # Efficiency = mapping_rate / time (higher is better)
        r['efficiency'] = r['mapping_rate'] / r['total_time_ms']
    
    best_efficiency = max(results, key=lambda x: x['efficiency'])
    
    print("⚡ BEST EFFICIENCY (mapping/time):")
    print(f"   Seeds: {best_efficiency['num_seeds']}, K: {best_efficiency['K']}")
    print(f"   Mapping: {best_efficiency['mapping_rate']:.1f}%")
    print(f"   Time: {best_efficiency['total_time_ms']:.1f}ms/read")
    print(f"   Efficiency: {best_efficiency['efficiency']:.3f} %/ms")
    print()
    
    # Tier recommendations
    print("="*80)
    print("ADAPTIVE TIER RECOMMENDATIONS")
    print("="*80 + "\n")
    
    # Find good tier candidates
    tier1_candidates = [r for r in results if 
                        r['mapping_rate'] >= 85 and 
                        r['total_time_ms'] <= 10]
    
    tier2_candidates = [r for r in results if 
                        r['mapping_rate'] >= 95]
    
    if tier1_candidates:
        tier1 = min(tier1_candidates, key=lambda x: x['total_time_ms'])
        print(f"🥇 TIER 1 (Fast - handles {tier1['mapping_rate']:.0f}% of reads):")
        print(f"   Seeds: {tier1['num_seeds']}, K: {tier1['K']}")
        print(f"   Time: {tier1['total_time_ms']:.1f}ms/read")
        print(f"   Exit criteria: chain_score ≥ {int(tier1['num_seeds'] * 0.5)}")
        print()
    
    if tier2_candidates:
        tier2 = min(tier2_candidates, key=lambda x: x['total_time_ms'])
        print(f"🥈 TIER 2 (Rescue - handles remaining {100 - tier1['mapping_rate']:.0f}%):")
        print(f"   Seeds: {tier2['num_seeds']}, K: {tier2['K']}")
        print(f"   Time: {tier2['total_time_ms']:.1f}ms/read")
        print(f"   Expected total: {tier2['mapping_rate']:.1f}% mapping")
        print()
        
        # Calculate average time for 2-tier
        if tier1_candidates:
            avg_time = (tier1['mapping_rate'] / 100 * tier1['total_time_ms'] + 
                       (100 - tier1['mapping_rate']) / 100 * tier2['total_time_ms'])
            print(f"   Average time (2-tier): {avg_time:.1f}ms/read")
            print(f"   Speedup vs Tier 2 alone: {tier2['total_time_ms'] / avg_time:.2f}x")
    
    print()
    
    # Quality metrics analysis
    print("="*80)
    print("QUALITY METRICS ANALYSIS")
    print("="*80 + "\n")
    
    # Correlation: anchors per seed vs mapping rate
    print("📊 Anchor Efficiency vs Mapping Rate:")
    for seeds in sorted(set(r['num_seeds'] for r in results)):
        seed_results = [r for r in results if r['num_seeds'] == seeds]
        print(f"\n   Seeds={seeds}:")
        for r in sorted(seed_results, key=lambda x: x['K']):
            efficiency = r['avg_anchor_efficiency']
            print(f"     K={r['K']}: {r['avg_anchors_per_seed']:.1f} anchors/seed "
                  f"({efficiency:.1%} of expected) → {r['mapping_rate']:.1f}% mapping")
    
    print("\n📊 Chain Quality Indicators:")
    for r in sorted(results, key=lambda x: -x['mapping_rate'])[:5]:
        print(f"   seeds={r['num_seeds']}, K={r['K']}: "
              f"score={r['avg_best_chain_score']:.1f}, "
              f"ambiguity={r['avg_chain_ambiguity']:.3f}, "
              f"chr_conc={r['avg_chr_concentration']:.3f}")
    
    print()


def main():
    """Analyze speed/accuracy sweep results"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze speed vs accuracy results')
    parser.add_argument('--results', default='results/speed_accuracy_sweep.csv',
                        help='Results CSV file')
    
    args = parser.parse_args()
    
    results_path = Path(__file__).parent / args.results
    
    if not results_path.exists():
        print(f"❌ Results file not found: {results_path}")
        print("Run speed_accuracy_sweep.py first!")
        return
    
    results = load_csv(results_path)
    analyze_speed_accuracy(results)


if __name__ == '__main__':
    main()
