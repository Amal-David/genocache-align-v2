#!/usr/bin/env python3
"""
Analyze Parameter Sweep Results

Visualizes and summarizes parameter sweep data to inform
intelligent tier design decisions.
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
                if key in ['num_seeds', 'K', 'tolerance', 'mapped', 'total_reads']:
                    row[key] = int(row[key])
                elif key not in ['similarity_threshold'] and row[key] != 'None':
                    try:
                        row[key] = float(row[key])
                    except:
                        pass
            results.append(row)
    return results


def analyze_k_sweep(results: List[Dict]):
    """Analyze K sweep results"""
    print("\n" + "="*80)
    print("K (NEIGHBORS) SWEEP ANALYSIS")
    print("="*80 + "\n")
    
    print(f"{'K':<8} {'Mapping %':<12} {'Time (ms)':<12} {'Anchors/Seed':<15} {'Chr Conc':<12}")
    print("-" * 80)
    
    for r in results:
        print(f"{r['K']:<8} {r['mapping_rate']:<12.1f} {r['avg_time_ms']:<12.1f} "
              f"{r['avg_anchor_density']:<15.2f} {r['avg_chr_concentration']:<12.3f}")
    
    # Find best
    best = max(results, key=lambda x: x['mapping_rate'])
    fastest = min(results, key=lambda x: x['avg_time_ms'])
    
    print("\n📊 Key Findings:")
    print(f"  • Best accuracy: K={best['K']} ({best['mapping_rate']:.1f}%)")
    print(f"  • Fastest: K={fastest['K']} ({fastest['avg_time_ms']:.1f}ms)")
    print(f"  • Anchor density scales with K: {results[0]['avg_anchor_density']:.1f} → {results[-1]['avg_anchor_density']:.1f}")
    
    # Recommendation
    if best['K'] == fastest['K']:
        print(f"\n💡 Recommendation: K={best['K']} (optimal for both speed and accuracy)")
    else:
        # Find sweet spot
        for r in results:
            if r['mapping_rate'] >= best['mapping_rate'] - 2:  # Within 2% of best
                print(f"\n💡 Recommendation: K={r['K']} (good accuracy {r['mapping_rate']:.1f}%, fast {r['avg_time_ms']:.1f}ms)")
                break


def analyze_tolerance_sweep(results: List[Dict]):
    """Analyze tolerance sweep results"""
    print("\n" + "="*80)
    print("TOLERANCE (C) SWEEP ANALYSIS")
    print("="*80 + "\n")
    
    print(f"{'C':<8} {'Mapping %':<12} {'Time (ms)':<12} {'Chains':<10} {'Ambiguity':<12}")
    print("-" * 80)
    
    for r in results:
        print(f"{r['tolerance']:<8} {r['mapping_rate']:<12.1f} {r['avg_time_ms']:<12.1f} "
              f"{r['avg_num_chains']:<10.2f} {r['avg_chain_ambiguity']:<12.3f}")
    
    # Find best
    best = max(results, key=lambda x: x['mapping_rate'])
    
    print("\n📊 Key Findings:")
    print(f"  • Best accuracy: C={best['tolerance']} ({best['mapping_rate']:.1f}%)")
    print(f"  • Chain count varies: {results[0]['avg_num_chains']:.1f} → {results[-1]['avg_num_chains']:.1f}")
    print(f"  • Higher C = more chains (more relaxed merging)")
    
    # Recommendation
    print(f"\n💡 Recommendation: C={best['tolerance']} for standard, "
          f"C={results[-1]['tolerance']} for difficult reads")


def analyze_similarity_sweep(results: List[Dict]):
    """Analyze similarity threshold sweep results"""
    print("\n" + "="*80)
    print("SIMILARITY THRESHOLD SWEEP ANALYSIS")
    print("="*80 + "\n")
    
    print(f"{'Threshold':<12} {'Mapping %':<12} {'Time (ms)':<12} {'Anchors/Seed':<15} {'Sim Mean':<12}")
    print("-" * 80)
    
    for r in results:
        threshold = r['similarity_threshold']
        print(f"{threshold:<12} {r['mapping_rate']:<12.1f} {r['avg_time_ms']:<12.1f} "
              f"{r['avg_anchor_density']:<15.2f} {r['avg_similarity_mean']:<12.3f}")
    
    # Find best
    best = max(results, key=lambda x: x['mapping_rate'])
    
    print("\n📊 Key Findings:")
    print(f"  • Best accuracy: threshold={best['similarity_threshold']} ({best['mapping_rate']:.1f}%)")
    
    # Check if filtering helps
    no_filter = [r for r in results if r['similarity_threshold'] == 'None'][0]
    with_filter = [r for r in results if r['similarity_threshold'] != 'None']
    
    if any(r['mapping_rate'] > no_filter['mapping_rate'] for r in with_filter):
        best_filtered = max(with_filter, key=lambda x: x['mapping_rate'])
        print(f"  • Filtering HELPS: {best_filtered['similarity_threshold']} gives "
              f"+{best_filtered['mapping_rate'] - no_filter['mapping_rate']:.1f}%")
    else:
        print(f"  • Filtering HURTS accuracy (loses {no_filter['mapping_rate'] - best['mapping_rate']:.1f}%)")
    
    # Recommendation
    if best['similarity_threshold'] == 'None':
        print(f"\n💡 Recommendation: No similarity filtering (keeps all anchors)")
    else:
        print(f"\n💡 Recommendation: Filter at {best['similarity_threshold']} "
              f"(+{best['mapping_rate'] - no_filter['mapping_rate']:.1f}% improvement)")


def generate_tier_recommendations(results_dir: Path):
    """Generate three-tier configuration recommendations"""
    print("\n" + "="*80)
    print("THREE-TIER CONFIGURATION RECOMMENDATIONS")
    print("="*80 + "\n")
    
    # Load all results
    k_file = results_dir / 'sweep_k.csv'
    c_file = results_dir / 'sweep_c.csv'
    sim_file = results_dir / 'sweep_similarity.csv'
    
    k_results = load_csv(k_file) if k_file.exists() else []
    c_results = load_csv(c_file) if c_file.exists() else []
    sim_results = load_csv(sim_file) if sim_file.exists() else []
    
    if not (k_results and c_results and sim_results):
        print("⚠️  Not all sweep results available yet")
        return
    
    # Find optimal values
    best_k = max(k_results, key=lambda x: x['mapping_rate'])['K']
    fast_k = min([r for r in k_results if r['mapping_rate'] >= 80], 
                 key=lambda x: x['avg_time_ms'])['K'] if any(r['mapping_rate'] >= 80 for r in k_results) else best_k
    
    best_c = max(c_results, key=lambda x: x['mapping_rate'])['tolerance']
    
    best_sim = max(sim_results, key=lambda x: x['mapping_rate'])
    use_sim_filter = best_sim['similarity_threshold'] != 'None'
    
    print("Based on parameter sweeps:\n")
    
    print("🚀 TIER 1 (Fast - 75-80% reads):")
    print(f"   Seeds: 16")
    print(f"   K: {fast_k} (optimized for speed)")
    print(f"   C: {best_c}")
    print(f"   Similarity: {best_sim['similarity_threshold'] if use_sim_filter else 'None'}")
    print(f"   Exit if: score ≥ 8 OR similarity > 0.80 OR ambiguity < 0.5")
    print(f"   Expected: ~250ms/read\n")
    
    print("⚡ TIER 2 (Standard - +15-17% reads):")
    print(f"   Seeds: 32")
    print(f"   K: {best_k}")
    print(f"   C: {int(best_c * 1.5)}")
    print(f"   Similarity: {0.65 if use_sim_filter else 'None'}")
    print(f"   Exit if: score ≥ 16 OR similarity > 0.70 OR ambiguity < 0.7")
    print(f"   Expected: ~500ms/read\n")
    
    print("🎯 TIER 3 (Deep - +2-4% reads):")
    print(f"   Seeds: 48")
    print(f"   K: {min(best_k + 8, 48)}")
    print(f"   C: {int(best_c * 2)}")
    print(f"   Similarity: {0.55 if use_sim_filter else 'None'}")
    print(f"   Exit: Never (final attempt)")
    print(f"   Expected: ~800ms/read\n")
    
    print("📊 Expected Performance:")
    print("   Total mapping: 97-99%")
    print("   Average time: ~350ms/read (vs 630ms for single-pass 32 seeds)")
    print("   Speedup: ~1.8x\n")


def main():
    """Analyze all sweep results"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze parameter sweep results')
    parser.add_argument('--results-dir', default='results', help='Results directory')
    
    args = parser.parse_args()
    
    results_dir = Path(__file__).parent / args.results_dir
    
    # Analyze each sweep
    k_file = results_dir / 'sweep_k.csv'
    if k_file.exists():
        k_results = load_csv(k_file)
        analyze_k_sweep(k_results)
    
    c_file = results_dir / 'sweep_c.csv'
    if c_file.exists():
        c_results = load_csv(c_file)
        analyze_tolerance_sweep(c_results)
    
    sim_file = results_dir / 'sweep_similarity.csv'
    if sim_file.exists():
        sim_results = load_csv(sim_file)
        analyze_similarity_sweep(sim_results)
    
    # Generate tier recommendations
    generate_tier_recommendations(results_dir)


if __name__ == '__main__':
    main()
