#!/usr/bin/env python3
"""
Parameter Sweep for NAL Pipeline

Tests different parameter combinations to understand:
1. K (neighbors per seed): Does more/fewer help?
2. C (tolerance): Optimal chaining tolerance?
3. Similarity threshold: Does filtering help?

Collects quality metrics for informed tier design.
"""

import sys
import time
import csv
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from seeding_nal import NALSeeding
from chaining_nal import NALChaining

class ParameterSweep:
    """
    Systematic parameter testing for NAL pipeline
    """
    
    def __init__(
        self,
        model_path: str,
        index_path: str,
        positions_path: str,
        device: str = 'cuda'
    ):
        """Initialize with NAL components"""
        print("Initializing Parameter Sweep...")
        
        # Initialize seeding (will modify K per test)
        self.model_path = model_path
        self.index_path = index_path
        self.positions_path = positions_path
        self.device = device
        
        self.seeder = None
        self.chainer = None
        
    def load_seeder(self, K: int):
        """Load seeder with specific K"""
        if self.seeder is not None:
            del self.seeder
        
        self.seeder = NALSeeding(
            model_path=self.model_path,
            index_path=self.index_path,
            positions_path=self.positions_path,
            seed_len=512,
            K=K,
            device=self.device
        )
    
    def load_chainer(self, tolerance: int):
        """Load chainer with specific tolerance"""
        self.chainer = NALChaining(
            tolerance=tolerance,
            min_chain_score=3,
            top_k=5
        )
    
    def compute_quality_metrics(
        self,
        anchors: List[Dict],
        chains: List[Dict],
        num_seeds: int
    ) -> Dict:
        """
        Compute anchor and chain quality metrics
        
        Metrics:
        - anchor_density: Anchors per seed
        - anchor_scatter: Std dev of positions per seed
        - similarity_mean: Average similarity score
        - similarity_min: Minimum similarity
        - chr_concentration: % anchors on top chromosome
        - chain_ambiguity: Score ratio (2nd / 1st)
        """
        metrics = {}
        
        if not anchors:
            return {
                'anchor_density': 0,
                'anchor_scatter': 0,
                'similarity_mean': 0,
                'similarity_min': 0,
                'chr_concentration': 0,
                'chain_ambiguity': 1.0,
                'num_chains': 0
            }
        
        # Anchor density
        metrics['anchor_density'] = len(anchors) / max(1, num_seeds)
        
        # Anchor scatter (std dev of ref positions)
        ref_positions = [a['ref_pos'] for a in anchors]
        metrics['anchor_scatter'] = float(np.std(ref_positions)) if len(ref_positions) > 1 else 0
        
        # Similarity distribution
        similarities = [a.get('similarity', 0) for a in anchors]
        metrics['similarity_mean'] = float(np.mean(similarities)) if similarities else 0
        metrics['similarity_min'] = float(np.min(similarities)) if similarities else 0
        
        # Chromosome concentration
        chr_counts = defaultdict(int)
        for a in anchors:
            chr_counts[a['ref_chr']] += 1
        if chr_counts:
            metrics['chr_concentration'] = max(chr_counts.values()) / len(anchors)
        else:
            metrics['chr_concentration'] = 0
        
        # Chain metrics
        metrics['num_chains'] = len(chains)
        if len(chains) >= 2:
            metrics['chain_ambiguity'] = chains[1]['score'] / max(1, chains[0]['score'])
        elif len(chains) == 1:
            metrics['chain_ambiguity'] = 0.0  # Clear winner
        else:
            metrics['chain_ambiguity'] = 1.0  # No chains
        
        return metrics
    
    def filter_anchors_by_similarity(
        self,
        anchors: List[Dict],
        threshold: float
    ) -> List[Dict]:
        """Filter anchors below similarity threshold"""
        if threshold is None:
            return anchors
        return [a for a in anchors if a.get('similarity', 0) >= threshold]
    
    def test_configuration(
        self,
        reads: List[Tuple[str, str]],
        num_seeds: int,
        K: int,
        tolerance: int,
        similarity_threshold: float = None
    ) -> Dict:
        """
        Test a single parameter configuration
        
        Returns:
            Results dict with metrics
        """
        # Load components with specific parameters
        if self.seeder is None or self.seeder.K != K:
            self.load_seeder(K)
        if self.chainer is None or self.chainer.tolerance != tolerance:
            self.load_chainer(tolerance)
        
        # Results
        mapped = 0
        total_time = 0
        all_metrics = []
        
        for read_id, read_seq in reads:
            start = time.time()
            
            # Get anchors
            anchors = self.seeder.get_anchors(read_seq, num_seeds=num_seeds, K=K)
            
            # Filter by similarity if specified
            if similarity_threshold is not None:
                anchors = self.filter_anchors_by_similarity(anchors, similarity_threshold)
            
            # Chain
            read_len = len(read_seq)
            chains, needs_rescue = self.chainer.chain_with_rescue_check(
                anchors, read_len=read_len, num_seeds=num_seeds, seed_len=512
            )
            
            # Compute metrics
            metrics = self.compute_quality_metrics(anchors, chains, num_seeds)
            all_metrics.append(metrics)
            
            # Check if mapped
            if chains:
                mapped += 1
            
            elapsed = time.time() - start
            total_time += elapsed
        
        # Aggregate results
        results = {
            'num_seeds': num_seeds,
            'K': K,
            'tolerance': tolerance,
            'similarity_threshold': similarity_threshold if similarity_threshold is not None else 'None',
            'mapped': mapped,
            'total_reads': len(reads),
            'mapping_rate': mapped / len(reads) * 100,
            'avg_time_ms': total_time / len(reads) * 1000,
        }
        
        # Average metrics
        for key in all_metrics[0].keys():
            values = [m[key] for m in all_metrics]
            results[f'avg_{key}'] = np.mean(values)
            results[f'std_{key}'] = np.std(values)
        
        return results
    
    def sweep_K(
        self,
        reads: List[Tuple[str, str]],
        num_seeds: int = 16,
        K_values: List[int] = [16, 24, 32, 40, 48],
        tolerance: int = 1000,
        output_file: str = 'results/sweep_k.csv'
    ):
        """Sweep K (neighbors) parameter"""
        print(f"\n{'='*80}")
        print(f"SWEEP: K (neighbors) with seeds={num_seeds}, tolerance={tolerance}")
        print(f"{'='*80}\n")
        
        results = []
        for K in K_values:
            print(f"Testing K={K}...")
            result = self.test_configuration(reads, num_seeds, K, tolerance)
            results.append(result)
            print(f"  → Mapping: {result['mapping_rate']:.1f}%, Time: {result['avg_time_ms']:.1f}ms")
        
        # Save to CSV
        self._save_results(results, output_file)
        print(f"\n✅ Results saved to {output_file}")
        
        return results
    
    def sweep_tolerance(
        self,
        reads: List[Tuple[str, str]],
        num_seeds: int = 16,
        K: int = 32,
        tolerance_values: List[int] = [500, 1000, 1500, 2000, 3000],
        output_file: str = 'results/sweep_c.csv'
    ):
        """Sweep tolerance (C) parameter"""
        print(f"\n{'='*80}")
        print(f"SWEEP: Tolerance (C) with seeds={num_seeds}, K={K}")
        print(f"{'='*80}\n")
        
        results = []
        for tolerance in tolerance_values:
            print(f"Testing C={tolerance}...")
            result = self.test_configuration(reads, num_seeds, K, tolerance)
            results.append(result)
            print(f"  → Mapping: {result['mapping_rate']:.1f}%, Time: {result['avg_time_ms']:.1f}ms")
        
        # Save to CSV
        self._save_results(results, output_file)
        print(f"\n✅ Results saved to {output_file}")
        
        return results
    
    def sweep_similarity(
        self,
        reads: List[Tuple[str, str]],
        num_seeds: int = 16,
        K: int = 32,
        tolerance: int = 1000,
        similarity_values: List[float] = [None, 0.5, 0.6, 0.7, 0.8],
        output_file: str = 'results/sweep_similarity.csv'
    ):
        """Sweep similarity threshold parameter"""
        print(f"\n{'='*80}")
        print(f"SWEEP: Similarity threshold with seeds={num_seeds}, K={K}, C={tolerance}")
        print(f"{'='*80}\n")
        
        results = []
        for sim_threshold in similarity_values:
            label = "None" if sim_threshold is None else f"{sim_threshold:.1f}"
            print(f"Testing similarity_threshold={label}...")
            result = self.test_configuration(reads, num_seeds, K, tolerance, sim_threshold)
            results.append(result)
            print(f"  → Mapping: {result['mapping_rate']:.1f}%, Time: {result['avg_time_ms']:.1f}ms")
        
        # Save to CSV
        self._save_results(results, output_file)
        print(f"\n✅ Results saved to {output_file}")
        
        return results
    
    def _save_results(self, results: List[Dict], output_file: str):
        """Save results to CSV"""
        output_path = Path(__file__).parent / output_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if not results:
            return
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)


def load_reads(fastq_path: str, max_reads: int = 100) -> List[Tuple[str, str]]:
    """Load reads from FASTQ"""
    reads = []
    with open(fastq_path) as f:
        while len(reads) < max_reads:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            read_id = header[1:].split()[0]
            reads.append((read_id, seq))
    return reads


def main():
    """Run parameter sweeps"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Parameter sweep for NAL pipeline')
    parser.add_argument('--model', required=True, help='Model path')
    parser.add_argument('--index', required=True, help='Index path')
    parser.add_argument('--positions', required=True, help='Positions path')
    parser.add_argument('--reads', required=True, help='FASTQ reads')
    parser.add_argument('--max-reads', type=int, default=100, help='Max reads to test')
    parser.add_argument('--sweep', choices=['k', 'c', 'similarity', 'all'], default='all', help='Which sweep to run')
    parser.add_argument('--device', default='cuda', help='Device')
    
    args = parser.parse_args()
    
    # Load reads
    print(f"Loading reads from {args.reads}...")
    reads = load_reads(args.reads, args.max_reads)
    print(f"  Loaded {len(reads)} reads\n")
    
    # Initialize sweep
    sweep = ParameterSweep(
        model_path=args.model,
        index_path=args.index,
        positions_path=args.positions,
        device=args.device
    )
    
    # Run sweeps
    if args.sweep in ['k', 'all']:
        sweep.sweep_K(reads)
    
    if args.sweep in ['c', 'all']:
        sweep.sweep_tolerance(reads)
    
    if args.sweep in ['similarity', 'all']:
        sweep.sweep_similarity(reads)
    
    print("\n" + "="*80)
    print("✅ All sweeps complete!")
    print("="*80)


if __name__ == '__main__':
    main()
