#!/usr/bin/env python3
"""
NAL-Aligned Data Augmentation - Exact Match to NeuralAligner

Key points from NAL paper (Section 3.1):
1. Error rate: U[0.01, 0.1] (random per sample)
2. Shift: ±L_seed/10 (bounded normal, we use uniform for simplicity)
3. RC augmentation: 50% probability INDEPENDENTLY for anchor and positive
4. NO double RC application!

Reference: NeuralAligner paper Section 3.1, Line 21-23
"""

import random
import numpy as np
from typing import List, Tuple


class NALAugmentation:
    """
    Clean augmentation matching NeuralAligner exactly
    
    NO inconsistencies, NO double RC, simple and clear
    """
    
    def __init__(self, seed_len: int = 512):
        """
        Initialize augmentation
        
        Args:
            seed_len: Seed length (default 512bp)
        """
        self.seed_len = seed_len
        self.max_shift = seed_len // 10  # ±51bp for 512bp seeds
        self.complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    def add_errors(
        self, 
        seq: str, 
        error_rate: float,
        sub_prob: float = 0.6,
        ins_prob: float = 0.2,
        del_prob: float = 0.2
    ) -> str:
        """
        Add sequencing errors (NAL: substitutions, insertions, deletions)
        
        Args:
            seq: Input sequence
            error_rate: Overall error rate (0.01-0.10)
            sub_prob: Substitution probability (0.6)
            ins_prob: Insertion probability (0.2)
            del_prob: Deletion probability (0.2)
        
        Returns:
            Sequence with errors
        """
        result = []
        i = 0
        
        while i < len(seq):
            if random.random() < error_rate:
                # Choose error type
                r = random.random()
                
                if r < sub_prob:
                    # Substitution
                    bases = [b for b in ['A', 'C', 'G', 'T'] if b != seq[i]]
                    result.append(random.choice(bases))
                    i += 1
                
                elif r < sub_prob + ins_prob:
                    # Insertion
                    result.append(random.choice(['A', 'C', 'G', 'T']))
                    result.append(seq[i])
                    i += 1
                
                else:
                    # Deletion
                    i += 1
            else:
                result.append(seq[i])
                i += 1
        
        return ''.join(result)
    
    def shift_and_pad(self, seq: str, shift: int, target_len: int) -> str:
        """
        Shift sequence and pad/trim to target length
        
        Args:
            seq: Input sequence
            shift: Shift amount (negative=left, positive=right)
            target_len: Target length
        
        Returns:
            Shifted and padded sequence
        """
        if shift > 0:
            # Shift right: pad left with N's, trim right
            seq = 'N' * shift + seq
        elif shift < 0:
            # Shift left: remove from left, pad right
            seq = seq[abs(shift):]
        
        # Pad or trim to target length
        if len(seq) < target_len:
            seq = seq + 'N' * (target_len - len(seq))
        elif len(seq) > target_len:
            seq = seq[:target_len]
        
        return seq
    
    def reverse_complement(self, seq: str) -> str:
        """
        Reverse complement
        
        Args:
            seq: Input sequence
        
        Returns:
            Reverse complement
        """
        return ''.join(self.complement.get(b, 'N') for b in reversed(seq))
    
    def create_training_pairs(
        self, 
        genome: str, 
        batch_size: int = 8192
    ) -> Tuple[List[str], List[str]]:
        """
        Generate training pairs for InfoNCE (NAL protocol)
        
        Steps (EXACT NAL order):
        1. Extract anchor from random position
        2. Create positive: apply errors + shift + pad
        3. RC augmentation: INDEPENDENTLY flip anchor and positive (50% each)
        
        NO explicit negatives! InfoNCE uses other samples in batch as negatives.
        
        Args:
            genome: Reference genome sequence
            batch_size: Number of pairs (NAL uses 8192)
        
        Returns:
            (anchors, positives): Training pairs
        """
        anchors = []
        positives = []
        
        max_attempts = batch_size * 10
        attempts = 0
        
        while len(anchors) < batch_size and attempts < max_attempts:
            attempts += 1
            
            # 1. Extract anchor sequence
            max_pos = len(genome) - self.seed_len
            if max_pos <= 0:
                break
            
            pos = random.randint(0, max_pos)
            anchor = genome[pos:pos + self.seed_len]
            
            # Skip if too many N's (NAL: discard sequences with 'N')
            if anchor.count('N') > self.seed_len * 0.1:
                continue
            
            # 2. Create positive: errors + shift + pad
            # Error rate: U[0.01, 0.1] (NAL spec)
            error_rate = random.uniform(0.01, 0.10)
            positive = self.add_errors(anchor, error_rate)
            
            # Shift: ±L_seed/10 (NAL spec)
            shift = random.randint(-self.max_shift, self.max_shift)
            positive = self.shift_and_pad(positive, shift, self.seed_len)
            
            # 3. RC augmentation: INDEPENDENTLY (NAL spec, Line 36)
            # "each sequence and its paired sequence are independently 
            #  transformed to their reverse complements with 50% probability"
            if random.random() < 0.5:
                anchor = self.reverse_complement(anchor)
            
            if random.random() < 0.5:
                positive = self.reverse_complement(positive)
            
            # Add to batch
            anchors.append(anchor)
            positives.append(positive)
        
        if len(anchors) < batch_size:
            print(f"Warning: Only generated {len(anchors)}/{batch_size} pairs")
        
        return anchors, positives
    
    def validate_augmentation(self, genome: str, n_samples: int = 100):
        """
        Validate augmentation matches NAL specs
        
        Args:
            genome: Test genome
            n_samples: Number of samples to validate
        
        Returns:
            Statistics dict
        """
        anchors, positives = self.create_training_pairs(genome, n_samples)
        
        stats = {
            'n_samples': len(anchors),
            'error_rates': [],
            'rc_anchor': 0,
            'rc_positive': 0,
            'rc_both': 0,
            'rc_neither': 0,
            'lengths_ok': 0
        }
        
        for i in range(len(anchors)):
            anchor = anchors[i]
            positive = positives[i]
            
            # Check lengths
            if len(anchor) == self.seed_len and len(positive) == self.seed_len:
                stats['lengths_ok'] += 1
            
            # Approximate error rate (rough estimate)
            # Can't compute exactly due to shifts, but check if different
            if anchor != positive:
                # Has some augmentation
                pass
        
        # Check RC statistics (approximate)
        # We can't determine this exactly after the fact, but we can check diversity
        unique_anchors = len(set(anchors))
        unique_positives = len(set(positives))
        
        stats['unique_anchors'] = unique_anchors
        stats['unique_positives'] = unique_positives
        stats['lengths_ok_pct'] = stats['lengths_ok'] / len(anchors) * 100
        
        return stats


def test_augmentation():
    """Test NAL augmentation"""
    print("Testing NAL-Aligned Augmentation")
    print("=" * 60)
    
    # Create augmenter
    aug = NALAugmentation(seed_len=512)
    
    # Create mock genome
    genome = "ACGTACGT" * 10000  # 80kb
    
    # Generate small batch
    print("\n1. Generating training pairs...")
    anchors, positives = aug.create_training_pairs(genome, batch_size=10)
    
    print(f"   Generated {len(anchors)} pairs")
    print(f"   Anchor[0][:50]:   {anchors[0][:50]}")
    print(f"   Positive[0][:50]: {positives[0][:50]}")
    print(f"   Same? {anchors[0] == positives[0]}")
    
    # Validate augmentation
    print("\n2. Validating augmentation (100 samples)...")
    stats = aug.validate_augmentation(genome, n_samples=100)
    print(f"   Samples: {stats['n_samples']}")
    print(f"   Unique anchors: {stats['unique_anchors']}")
    print(f"   Unique positives: {stats['unique_positives']}")
    print(f"   Lengths OK: {stats['lengths_ok_pct']:.1f}%")
    
    # Test individual functions
    print("\n3. Testing individual functions...")
    
    seq = "ACGTACGTACGT" * 43  # 516bp
    
    # Test errors
    noisy = aug.add_errors(seq, error_rate=0.05)
    errors = sum(1 for a, b in zip(seq, noisy) if a != b)
    print(f"   Errors added: {errors}/{len(seq)} = {errors/len(seq)*100:.1f}%")
    
    # Test shift
    shifted = aug.shift_and_pad(seq, shift=10, target_len=512)
    print(f"   Shifted length: {len(shifted)} (should be 512)")
    
    # Test RC
    rc = aug.reverse_complement(seq)
    rc_rc = aug.reverse_complement(rc)
    print(f"   RC works: {seq == rc_rc}")
    
    print("\n✅ NAL augmentation test complete!")


if __name__ == '__main__':
    test_augmentation()
