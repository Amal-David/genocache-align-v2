#!/usr/bin/env python3
"""
NAL-Aligned Chaining - Exact NeuralAligner Protocol

Key points from NAL paper (Section 3.4):
- Relaxed chaining for inexact anchors
- Positional tolerance C (default 1000bp)
- Select diagonal stripe with slope ±1
- Score = number of anchors (not span length!)
- Return top-K chains

Reference: NeuralAligner paper Section 3.4, Equation 2
"""

import numpy as np
from typing import List, Dict, Tuple
from collections import defaultdict


class NALChaining:
    """
    NAL-aligned chaining following paper Section 3.4
    
    Handles inexact anchors with relaxed colinearity constraints
    """
    
    def __init__(
        self,
        tolerance: int = 1000,  # NAL: C parameter (positional tolerance)
        min_chain_score: int = 3,  # Minimum anchors in chain
        top_k: int = 5  # Return top-K chains
    ):
        """
        Initialize NAL chaining
        
        Args:
            tolerance: Positional tolerance C (NAL Eq. 2) and C_b (Alg. 1)
                      Paper distinguishes between:
                      - C: colinearity tolerance in Eq. 2 (~1000bp)
                      - C_b: bias tolerance for merging in Alg. 1
                      We use same value for both (reasonable simplification)
            min_chain_score: Minimum chain score to keep (our design: 3)
            top_k: Number of top chains to return (NAL: top-K)
        """
        self.tolerance = tolerance
        self.min_chain_score = min_chain_score
        self.top_k = top_k
    
    def chain_anchors(
        self,
        anchors: List[Dict],
        read_len: int,
        seed_len: int = 512
    ) -> List[Dict]:
        """
        Chain anchors using NAL algorithm (Section 3.4, Equation 2)
        
        NAL chaining:
        - For each strand (+/-), find diagonal stripe
        - Stripe has slope 1 (positive) or -1 (negative)
        - Width = 2×tolerance
        - Score = number of anchors in stripe
        - Return top-K chains
        
        Args:
            anchors: List of anchor dicts
            read_len: Read length
            seed_len: Seed length
        
        Returns:
            List of chain dicts with:
                - ref_chr: Chromosome
                - ref_pos: Leftmost genomic position
                - strand: '+' or '-'
                - score: Number of anchors in chain
                - anchors: List of anchors in chain
        """
        if not anchors:
            return []
        
        # Group anchors by chromosome
        chr_anchors = defaultdict(list)
        for anchor in anchors:
            chr_anchors[anchor['ref_chr']].append(anchor)
        
        all_chains = []
        
        # Process each chromosome
        for ref_chr, chr_anchor_list in chr_anchors.items():
            # Try both strands
            for strand in ['+', '-']:
                chains = self._find_chains_on_strand(
                    chr_anchor_list,
                    ref_chr,
                    strand,
                    read_len,
                    seed_len
                )
                all_chains.extend(chains)
        
        # Sort by score (descending)
        all_chains.sort(key=lambda x: x['score'], reverse=True)
        
        # Return top-K chains
        return all_chains[:self.top_k]
    
    def _find_chains_on_strand(
        self,
        anchors: List[Dict],
        ref_chr: str,
        strand: str,
        read_len: int,
        seed_len: int
    ) -> List[Dict]:
        """
        Find chains on specific strand using NAL Algorithm 1
        
        Follows paper Section A.4, Algorithm 1 (Vectorized chaining):
        1. Adjust Y for strand (Y = Y - X for +, Y = Y + X + L_seed - L_read for -)
        2. Sort anchors by adjusted Y
        3. Merge neighbors within C_b (bias tolerance)
        4. Use similarity D to pick best anchor in merge
        5. Score = count of anchors in merged stripe
        
        Args:
            anchors: Anchors on this chromosome
            ref_chr: Chromosome name
            strand: '+' or '-'
            read_len: Read length
            seed_len: Seed length
        
        Returns:
            List of chains
        """
        if not anchors:
            return []
        
        # Step 1: Adjust Y for strand (Algorithm 1, line 3)
        adjusted_anchors = []
        for anchor in anchors:
            y_ij = anchor['ref_pos']  # Genomic position
            x_i = anchor['seed_pos']  # Position in read
            similarity = anchor.get('similarity', 0.0)
            
            if strand == '+':
                # Positive strand: Y = Y - X
                adjusted_y = y_ij - x_i
            else:
                # Negative strand: Y = Y + X + L_seed - L_read
                adjusted_y = y_ij + x_i + seed_len - read_len
            
            adjusted_anchors.append({
                'original': anchor,
                'adjusted_y': adjusted_y,
                'similarity': similarity,
                'score': 1  # Initial score (Algorithm 1, line 1)
            })
        
        # Step 2: Sort by adjusted Y (Algorithm 1, line 5)
        adjusted_anchors.sort(key=lambda a: a['adjusted_y'])
        
        # Step 3: Merge neighbors within C_b (Algorithm 1, lines 7-12)
        # C_b is the bias tolerance for merging
        C_b = self.tolerance  # Using same tolerance as bias tolerance
        
        # Track which anchors belong to which stripe
        # Each stripe leader keeps list of all merged anchors
        for anchor in adjusted_anchors:
            anchor['merged_anchors'] = [anchor['original']]
        
        # Mark anchors to keep (not merged away)
        active = [True] * len(adjusted_anchors)
        
        for i in range(len(adjusted_anchors)):
            if not active[i]:
                continue
            
            # Merge all neighbors within C_b into i
            y_i = adjusted_anchors[i]['adjusted_y']
            
            for j in range(i + 1, len(adjusted_anchors)):
                if not active[j]:
                    continue
                
                y_j = adjusted_anchors[j]['adjusted_y']
                
                # Stop if beyond C_b range (sorted, so no more matches possible)
                if y_j - y_i > C_b:
                    break
                
                # Merge j into i (Algorithm 1, lines 8-11)
                # Pick max similarity as stripe axis
                if adjusted_anchors[j]['similarity'] > adjusted_anchors[i]['similarity']:
                    # j has better similarity, use it as axis
                    adjusted_anchors[i]['adjusted_y'] = y_j
                    adjusted_anchors[i]['similarity'] = adjusted_anchors[j]['similarity']
                
                # Accumulate score (count)
                adjusted_anchors[i]['score'] += adjusted_anchors[j]['score']
                
                # Merge anchor lists
                adjusted_anchors[i]['merged_anchors'].extend(adjusted_anchors[j]['merged_anchors'])
                
                # Mark j as merged (inactive)
                active[j] = False
        
        # Step 4: Collect active chains (Algorithm 1, line 13)
        chains = []
        for i, anchor_data in enumerate(adjusted_anchors):
            if not active[i]:
                continue
            
            score = anchor_data['score']
            
            if score >= self.min_chain_score:
                # Collect all anchors in this stripe
                chain_anchors = anchor_data['merged_anchors']
                
                chains.append({
                    'ref_chr': ref_chr,
                    'ref_pos': int(anchor_data['adjusted_y']),
                    'strand': strand,
                    'score': score,
                    'anchors': chain_anchors,
                    'num_seeds': len(set(a['seed_idx'] for a in chain_anchors)),
                    'best_similarity': anchor_data['similarity']
                })
        
        # Step 5: Sort by score descending (Algorithm 1, line 13)
        chains.sort(key=lambda c: c['score'], reverse=True)
        
        return chains
    
    def chain_with_rescue_check(
        self,
        anchors: List[Dict],
        read_len: int,
        num_seeds: int,
        seed_len: int = 512
    ) -> Tuple[List[Dict], bool]:
        """
        Chain anchors and determine if rescue is needed (NAL Section 3.5)
        
        NAL rescue strategy (NAL A.5: accuracy mode):
        - First iteration: use 7 seeds (224 anchors with K=32)
        - Check if chain score < threshold (e.g., num_seeds/2)
        - If low score or ambiguous top chains → rescue needed
        - Rescue: add more seeds (13 total, 416 anchors with K=32)
        
        Args:
            anchors: List of anchors
            read_len: Read length
            num_seeds: Number of seeds used
            seed_len: Seed length
        
        Returns:
            (chains, needs_rescue): chains and rescue flag
        """
        chains = self.chain_anchors(anchors, read_len, seed_len)
        
        if not chains:
            return [], True  # No chains → rescue
        
        best_chain = chains[0]
        best_score = best_chain['score']
        
        # Check 1: Low score (NAL: chain contains few anchors)
        # Paper: "chains that fall below a score threshold"
        # Our design: threshold = num_seeds/2 (reasonable heuristic)
        score_threshold = num_seeds // 2
        if best_score < score_threshold:
            return chains, True
        
        # Check 2: Ambiguous top chains (NAL: similar scores)
        # Paper: "show ambiguous top-K results"
        # Our design: second score ≥ 80% of best (reasonable threshold)
        if len(chains) > 1:
            second_score = chains[1]['score']
            # If top 2 chains have similar scores → ambiguous
            if second_score >= best_score * 0.8:  # Within 20%
                return chains, True
        
        # Check 3: Maximum score achieved (NAL: explicit criterion)
        # Paper: "best chain achieves maximum score, equals number of seeds"
        if best_score == num_seeds:
            return chains, False  # Perfect chain, no rescue
        
        # Default: moderate score, no rescue needed
        return chains, False


def test_chaining():
    """Test NAL chaining"""
    print("="* 80)
    print("NAL Chaining Test")
    print("=" * 80)
    
    # Create mock anchors
    anchors = []
    
    # Simulate a good colinear chain on positive strand
    # Read: 1000bp, seeds at positions 0, 250, 500, 750
    # Ref: chr1, starting at position 1000000
    for i in range(4):
        seed_pos = i * 250
        ref_pos = 1000000 + seed_pos + np.random.randint(-50, 50)  # Small noise
        
        for j in range(8):  # 8 anchors per seed (K=32 total)
            anchors.append({
                'seed_idx': i,
                'seed_pos': seed_pos,
                'ref_chr': 'chr1',
                'ref_pos': ref_pos + j * 100,  # Some spread
                'similarity': 0.9 + np.random.rand() * 0.1,
                'anchor_idx': j
            })
    
    # Add some noise anchors (wrong positions)
    for i in range(5):
        anchors.append({
            'seed_idx': 0,
            'seed_pos': 0,
            'ref_chr': 'chr1',
            'ref_pos': 5000000 + i * 1000,  # Far away
            'similarity': 0.7,
            'anchor_idx': i
        })
    
    print(f"\nTest data:")
    print(f"  Total anchors: {len(anchors)}")
    print(f"  Read length: 1000bp")
    print(f"  Seeds: 4")
    
    # Chain
    chainer = NALChaining(tolerance=1000, min_chain_score=3, top_k=5)
    
    chains = chainer.chain_anchors(anchors, read_len=1000, seed_len=512)
    
    print(f"\nResults:")
    print(f"  Chains found: {len(chains)}")
    
    for i, chain in enumerate(chains):
        print(f"\nChain {i+1}:")
        print(f"  Chr: {chain['ref_chr']}")
        print(f"  Pos: {chain['ref_pos']:,}")
        print(f"  Strand: {chain['strand']}")
        print(f"  Score: {chain['score']} anchors")
        print(f"  Seeds: {chain['num_seeds']} unique")
    
    # Test rescue check
    print(f"\n" + "-" * 80)
    print("Testing rescue strategy:")
    
    for num_seeds in [4, 5]:
        chains, needs_rescue = chainer.chain_with_rescue_check(
            anchors, read_len=1000, num_seeds=num_seeds, seed_len=512
        )
        
        print(f"\n  With {num_seeds} seeds:")
        print(f"    Best score: {chains[0]['score'] if chains else 0}")
        print(f"    Needs rescue: {needs_rescue}")
    
    print("\n✅ Chaining test complete!")


if __name__ == '__main__':
    test_chaining()
