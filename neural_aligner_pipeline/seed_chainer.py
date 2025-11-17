#!/usr/bin/env python3
"""
Phase 2: Seed chaining using dynamic programming
Links co-linear seeds to improve accuracy 73.5% → 85-90%
"""

import numpy as np

class SeedChainer:
    def __init__(self, max_gap=10000, max_deviation=1000):
        """
        Initialize seed chainer
        
        Args:
            max_gap: Maximum distance between consecutive seeds (bp)
            max_deviation: Maximum deviation from co-linearity (bp)
        """
        self.max_gap = max_gap
        self.max_deviation = max_deviation
    
    def chain_seeds(self, seeds):
        """
        Chain seeds using dynamic programming
        
        Args:
            seeds: List of seed hit dictionaries with:
                - seed_read_pos: Position in read
                - seed_ref_pos: Position in reference
                - score: Similarity score
        
        Returns:
            Best chain of co-linear seeds
        """
        if not seeds:
            return []
        
        # Sort seeds by position in read
        seeds = sorted(seeds, key=lambda s: s['seed_read_pos'])
        
        n = len(seeds)
        
        # DP arrays
        best_score = [s['score'] for s in seeds]  # Best score ending at seed i
        predecessor = [-1] * n  # Previous seed in best chain
        
        # Dynamic programming: for each seed, find best predecessor
        for i in range(1, n):
            for j in range(i):
                # Check if seed j and seed i are co-linear
                if self.are_colinear(seeds[j], seeds[i]):
                    # Score = previous chain score + current seed score
                    score = best_score[j] + seeds[i]['score']
                    
                    if score > best_score[i]:
                        best_score[i] = score
                        predecessor[i] = j
        
        # Traceback to get best chain
        best_end = max(range(n), key=lambda i: best_score[i])
        
        chain = []
        idx = best_end
        
        while idx != -1:
            chain.append(seeds[idx])
            idx = predecessor[idx]
        
        chain.reverse()
        
        return chain
    
    def are_colinear(self, seed1, seed2):
        """
        Check if two seeds are co-linear (consistent alignment)
        
        Co-linear means:
        - Both seeds advance in same direction
        - Distance in read ≈ distance in reference (within tolerance)
        """
        # Distance in read
        read_dist = seed2['seed_read_pos'] - seed1['seed_read_pos']
        
        # Distance in reference
        ref_dist = seed2['seed_ref_pos'] - seed1['seed_ref_pos']
        
        # Both must be positive (advancing forward)
        if read_dist <= 0 or ref_dist <= 0:
            return False
        
        # Gap must not be too large
        if read_dist > self.max_gap or ref_dist > self.max_gap:
            return False
        
        # Deviation (accounts for indels)
        deviation = abs(read_dist - ref_dist)
        
        return deviation < self.max_deviation
    
    def chain_clusters(self, clusters, allow_single_seed=True):
        """
        Chain seeds within each cluster and return best alignment
        
        Args:
            clusters: List of cluster dictionaries from multi-seeder
            allow_single_seed: If True, return single-seed hits if no chain found
        
        Returns:
            Best cluster with chained seeds
        """
        best_cluster = None
        best_chain_score = 0
        best_single_seed = None
        best_single_score = 0
        
        for cluster in clusters:
            # Get all hits in this cluster
            seeds = cluster['hits']
            
            # Chain seeds
            chain = self.chain_seeds(seeds)
            
            if chain:
                # Calculate chain quality metrics
                chain_score = sum(s['score'] for s in chain)
                chain_length = len(chain)
                
                # Coverage: how much of read is covered by chain
                if chain_length > 1:
                    read_coverage = (chain[-1]['seed_read_pos'] - 
                                   chain[0]['seed_read_pos'])
                else:
                    read_coverage = 0
                
                # Combined score: favor longer chains with high scores
                combined_score = chain_score * (1 + 0.1 * chain_length)
                
                if combined_score > best_chain_score:
                    best_chain_score = combined_score
                    best_cluster = {
                        'ref_pos': chain[0]['seed_ref_pos'] - chain[0]['seed_read_pos'],
                        'chain': chain,
                        'chain_length': chain_length,
                        'chain_score': chain_score,
                        'read_coverage': read_coverage,
                        'avg_score': chain_score / chain_length if chain_length > 0 else 0
                    }
                
                # Track best single seed as fallback
                if chain_length == 1 and chain[0]['score'] > best_single_score:
                    best_single_score = chain[0]['score']
                    best_single_seed = {
                        'ref_pos': chain[0]['seed_ref_pos'] - chain[0]['seed_read_pos'],
                        'chain': chain,
                        'chain_length': 1,
                        'chain_score': chain[0]['score'],
                        'read_coverage': 0,
                        'avg_score': chain[0]['score']
                    }
        
        # Always return something - chaining is for refinement, not rejection
        if best_cluster and best_cluster['chain_length'] >= 2:
            # Have a good chain
            return best_cluster
        elif best_single_seed:
            # Fallback to best single seed
            return best_single_seed
        elif clusters:
            # Last resort: return top cluster even without good chain
            # Use the seed with highest score
            top_cluster = clusters[0]
            if top_cluster['hits']:
                best_hit = max(top_cluster['hits'], key=lambda h: h['score'])
                return {
                    'ref_pos': best_hit['seed_ref_pos'] - best_hit['seed_read_pos'],
                    'chain': [best_hit],
                    'chain_length': 1,
                    'chain_score': best_hit['score'],
                    'read_coverage': 0,
                    'avg_score': best_hit['score']
                }
        
        return None


class ChainedMultiSeedAligner:
    """Multi-seed aligner with chaining support"""
    
    def __init__(self, multi_seed_aligner, max_gap=10000, max_deviation=1000):
        """
        Initialize chained aligner
        
        Args:
            multi_seed_aligner: MultiSeedAligner instance
            max_gap: Maximum gap between seeds
            max_deviation: Maximum deviation for co-linearity
        """
        self.aligner = multi_seed_aligner
        self.chainer = SeedChainer(max_gap=max_gap, max_deviation=max_deviation)
    
    def align_read(self, read_seq, n_seeds=5, top_k=50):
        """
        Align read using multi-seeding + chaining
        
        Args:
            read_seq: Read sequence
            n_seeds: Number of seeds to extract
            top_k: Top K candidates per seed
        
        Returns:
            Best chained alignment
        """
        # Get clusters from multi-seeder
        clusters = self.aligner.align_read(
            read_seq, 
            n_seeds=n_seeds, 
            top_k=top_k,
            min_score=0.5,
            cluster_dist=500
        )
        
        if not clusters:
            return None
        
        # Chain seeds within clusters
        best_alignment = self.chainer.chain_clusters(clusters)
        
        return best_alignment


def main():
    """Test seed chaining"""
    print("="*70)
    print("Seed Chaining - Phase 2")
    print("="*70)
    
    # Test data: synthetic seeds that should chain
    test_seeds = [
        {'seed_read_pos': 0, 'seed_ref_pos': 1000, 'score': 0.9},
        {'seed_read_pos': 500, 'seed_ref_pos': 1500, 'score': 0.85},
        {'seed_read_pos': 1000, 'seed_ref_pos': 2000, 'score': 0.88},
        {'seed_read_pos': 1500, 'seed_ref_pos': 2500, 'score': 0.92},
        # Add some noise
        {'seed_read_pos': 200, 'seed_ref_pos': 5000, 'score': 0.7},  # Off-target
        {'seed_read_pos': 800, 'seed_ref_pos': 10000, 'score': 0.65},  # Off-target
    ]
    
    chainer = SeedChainer(max_gap=10000, max_deviation=1000)
    
    print(f"\nTest seeds: {len(test_seeds)}")
    for i, s in enumerate(test_seeds):
        print(f"  {i+1}. Read pos: {s['seed_read_pos']}, Ref pos: {s['seed_ref_pos']}, Score: {s['score']:.2f}")
    
    # Chain seeds
    chain = chainer.chain_seeds(test_seeds)
    
    print(f"\nChained seeds: {len(chain)}")
    for i, s in enumerate(chain):
        print(f"  {i+1}. Read pos: {s['seed_read_pos']}, Ref pos: {s['seed_ref_pos']}, Score: {s['score']:.2f}")
    
    total_score = sum(s['score'] for s in chain)
    print(f"\nChain score: {total_score:.2f}")
    print(f"Chain length: {len(chain)}")
    
    # Check co-linearity
    print(f"\nCo-linearity check:")
    for i in range(len(chain) - 1):
        s1, s2 = chain[i], chain[i+1]
        read_dist = s2['seed_read_pos'] - s1['seed_read_pos']
        ref_dist = s2['seed_ref_pos'] - s1['seed_ref_pos']
        deviation = abs(read_dist - ref_dist)
        print(f"  {i+1} → {i+2}: read_dist={read_dist}, ref_dist={ref_dist}, deviation={deviation}")
    
    print("\n" + "="*70)
    print("✓ Chaining test complete")
    print("="*70)


if __name__ == "__main__":
    main()
