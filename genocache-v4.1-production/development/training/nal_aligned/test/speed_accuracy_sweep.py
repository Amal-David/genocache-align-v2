#!/usr/bin/env python3
"""
Speed vs Accuracy Sweep - Finding the Sweet Spot

Tests combinations of seeds and K to find optimal speed/accuracy trade-off:
- Seeds: [2, 6, 8, 12, 16] (low to moderate coverage)
- K: [32, 48, 64] (standard to high neighbors)
- C: 1000 (fixed, proven optimal)

Tracks detailed timing and quality metrics for informed decisions.
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


class DetailedParameterSweep:
    """
    Comprehensive parameter testing with timing breakdown
    """
    
    def __init__(
        self,
        model_path: str,
        index_path: str,
        positions_path: str,
        device: str = 'cuda'
    ):
        """Initialize with NAL components"""
        print("Initializing Speed/Accuracy Sweep...")
        
        self.model_path = model_path
        self.index_path = index_path
        self.positions_path = positions_path
        self.device = device
        
        self.seeder = None
        self.chainer = None
        self.current_K = None
        
    def load_seeder(self, K: int):
        """Load seeder with specific K"""
        if self.current_K == K and self.seeder is not None:
            return  # Already loaded
        
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
        self.current_K = K
    
    def load_chainer(self, tolerance: int):
        """Load chainer with specific tolerance"""
        if self.chainer is None or self.chainer.tolerance != tolerance:
            self.chainer = NALChaining(
                tolerance=tolerance,
                min_chain_score=3,
                top_k=5
            )
    
    def compute_detailed_metrics(
        self,
        anchors: List[Dict],
        chains: List[Dict],
        num_seeds: int,
        K: int
    ) -> Dict:
        """
        Compute comprehensive quality metrics
        
        Key metrics for decision making:
        - Coverage: anchors per seed (expect ~K)
        - Quality: similarity scores
        - Specificity: chromosome concentration
        - Confidence: chain ambiguity
        - Chain strength: score relative to seeds
        """
        metrics = {}
        
        if not anchors:
            return {
                'total_anchors': 0,
                'anchors_per_seed': 0.0,
                'expected_anchors': num_seeds * K,
                'anchor_efficiency': 0.0,
                'similarity_mean': 0.0,
                'similarity_min': 0.0,
                'similarity_std': 0.0,
                'chr_concentration': 0.0,
                'chr_diversity': 0,
                'anchor_scatter': 0.0,
                'num_chains': 0,
                'best_chain_score': 0,
                'chain_score_ratio': 0.0,
                'chain_ambiguity': 1.0,
            }
        
        # Anchor coverage
        metrics['total_anchors'] = len(anchors)
        metrics['anchors_per_seed'] = len(anchors) / max(1, num_seeds)
        metrics['expected_anchors'] = num_seeds * K
        metrics['anchor_efficiency'] = len(anchors) / max(1, num_seeds * K)
        
        # Similarity distribution
        similarities = [a.get('similarity', 0) for a in anchors]
        metrics['similarity_mean'] = float(np.mean(similarities)) if similarities else 0
        metrics['similarity_min'] = float(np.min(similarities)) if similarities else 0
        metrics['similarity_std'] = float(np.std(similarities)) if similarities else 0
        
        # Chromosome specificity
        chr_counts = defaultdict(int)
        for a in anchors:
            chr_counts[a['ref_chr']] += 1
        
        if chr_counts:
            metrics['chr_concentration'] = max(chr_counts.values()) / len(anchors)
            metrics['chr_diversity'] = len(chr_counts)
        else:
            metrics['chr_concentration'] = 0
            metrics['chr_diversity'] = 0
        
        # Anchor scatter (positional spread)
        ref_positions = [a['ref_pos'] for a in anchors]
        metrics['anchor_scatter'] = float(np.std(ref_positions)) if len(ref_positions) > 1 else 0
        
        # Chain metrics
        metrics['num_chains'] = len(chains)
        
        if chains:
            metrics['best_chain_score'] = chains[0]['score']
            metrics['chain_score_ratio'] = chains[0]['score'] / max(1, num_seeds)
            
            if len(chains) >= 2:
                metrics['chain_ambiguity'] = chains[1]['score'] / max(1, chains[0]['score'])
            else:
                metrics['chain_ambiguity'] = 0.0  # Clear winner
        else:
            metrics['best_chain_score'] = 0
            metrics['chain_score_ratio'] = 0.0
            metrics['chain_ambiguity'] = 1.0  # No clear answer
        
        return metrics
    
    def test_configuration_detailed(
        self,
        reads: List[Tuple[str, str]],
        num_seeds: int,
        K: int,
        tolerance: int = 1000
    ) -> Dict:
        """
        Test configuration with detailed timing breakdown
        
        Returns:
            Results dict with metrics and timing
        """
        # Load components
        self.load_seeder(K)
        self.load_chainer(tolerance)
        
        # Results
        mapped = 0
        time_encoding = 0
        time_faiss = 0
        time_chaining = 0
        all_metrics = []
        
        for read_id, read_seq in reads:
            # Timing: Encoding
            t0 = time.time()
            anchors = self.seeder.get_anchors(read_seq, num_seeds=num_seeds, K=K)
            t1 = time.time()
            
            # Note: FAISS query time is included in get_anchors
            # We can approximate: encoding vs FAISS
            time_encoding += (t1 - t0) * 0.7  # ~70% encoding
            time_faiss += (t1 - t0) * 0.3     # ~30% FAISS
            
            # Timing: Chaining
            t2 = time.time()
            read_len = len(read_seq)
            chains, needs_rescue = self.chainer.chain_with_rescue_check(
                anchors, read_len=read_len, num_seeds=num_seeds, seed_len=512
            )
            t3 = time.time()
            time_chaining += (t3 - t2)
            
            # Compute metrics
            metrics = self.compute_detailed_metrics(anchors, chains, num_seeds, K)
            all_metrics.append(metrics)
            
            # Check if mapped
            if chains:
                mapped += 1
        
        # Aggregate results
        results = {
            'num_seeds': num_seeds,
            'K': K,
            'tolerance': tolerance,
            'mapped': mapped,
            'total_reads': len(reads),
            'mapping_rate': mapped / len(reads) * 100,
            'total_time_ms': (time_encoding + time_faiss + time_chaining) / len(reads) * 1000,
            'encoding_time_ms': time_encoding / len(reads) * 1000,
            'faiss_time_ms': time_faiss / len(reads) * 1000,
            'chaining_time_ms': time_chaining / len(reads) * 1000,
        }
        
        # Average metrics
        for key in all_metrics[0].keys():
            values = [m[key] for m in all_metrics]
            results[f'avg_{key}'] = np.mean(values)
            results[f'std_{key}'] = np.std(values)
        
        return results
    
    def run_full_sweep(
        self,
        reads: List[Tuple[str, str]],
        seed_values: List[int] = [2, 6, 8, 12, 16],
        K_values: List[int] = [32, 48, 64],
        tolerance: int = 1000,
        output_file: str = 'results/speed_accuracy_sweep.csv'
    ):
        """
        Run full speed/accuracy sweep
        
        Tests all combinations of seeds and K
        """
        print(f"\n{'='*80}")
        print(f"SPEED vs ACCURACY SWEEP")
        print(f"{'='*80}\n")
        print(f"Testing {len(seed_values)} seed counts × {len(K_values)} K values = {len(seed_values) * len(K_values)} configurations")
        print(f"Seeds: {seed_values}")
        print(f"K: {K_values}")
        print(f"Tolerance: {tolerance} (fixed)")
        print(f"Reads: {len(reads)}")
        print()
        
        results = []
        total_configs = len(seed_values) * len(K_values)
        current = 0
        
        for num_seeds in seed_values:
            for K in K_values:
                current += 1
                print(f"[{current}/{total_configs}] Testing seeds={num_seeds}, K={K}...")
                
                result = self.test_configuration_detailed(reads, num_seeds, K, tolerance)
                results.append(result)
                
                print(f"  → Mapping: {result['mapping_rate']:.1f}%, "
                      f"Time: {result['total_time_ms']:.1f}ms "
                      f"(enc: {result['encoding_time_ms']:.1f}ms, "
                      f"faiss: {result['faiss_time_ms']:.1f}ms, "
                      f"chain: {result['chaining_time_ms']:.1f}ms), "
                      f"Anchors: {result['avg_anchors_per_seed']:.1f}/seed")
        
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
    """Run speed/accuracy sweep"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Speed vs accuracy sweep for NAL pipeline')
    parser.add_argument('--model', required=True, help='Model path')
    parser.add_argument('--index', required=True, help='Index path')
    parser.add_argument('--positions', required=True, help='Positions path')
    parser.add_argument('--reads', required=True, help='FASTQ reads')
    parser.add_argument('--max-reads', type=int, default=100, help='Max reads to test')
    parser.add_argument('--seeds', type=int, nargs='+', default=[2, 6, 8, 12, 16], 
                        help='Seed counts to test')
    parser.add_argument('--K', type=int, nargs='+', default=[32, 48, 64],
                        help='K values to test')
    parser.add_argument('--tolerance', type=int, default=1000, help='Chaining tolerance')
    parser.add_argument('--device', default='cuda', help='Device')
    
    args = parser.parse_args()
    
    # Load reads
    print(f"Loading reads from {args.reads}...")
    reads = load_reads(args.reads, args.max_reads)
    print(f"  Loaded {len(reads)} reads\n")
    
    # Initialize sweep
    sweep = DetailedParameterSweep(
        model_path=args.model,
        index_path=args.index,
        positions_path=args.positions,
        device=args.device
    )
    
    # Run sweep
    sweep.run_full_sweep(
        reads=reads,
        seed_values=args.seeds,
        K_values=args.K,
        tolerance=args.tolerance
    )
    
    print("\n" + "="*80)
    print("✅ Speed/Accuracy sweep complete!")
    print("="*80)


if __name__ == '__main__':
    main()
