#!/usr/bin/env python3
"""
GenoCache V4 - Data Augmentation Pipeline

Implements NeuralAligner-proven augmentation strategies:
1. Error injection (1-10% - mimics sequencing errors)
2. Position shift (±50bp - translation continuity)
3. RC augmentation (50% - strand invariance)

Reference: NeuralAligner paper (ICLR 2026 submission)
"""

import random
import numpy as np
from typing import Tuple, List, Dict, Optional


class DataAugmentation:
    """
    NeuralAligner-proven augmentation strategies
    """
    
    def __init__(self, seed_len: int = 512):
        """
        Initialize augmentation pipeline
        
        Args:
            seed_len: Length of DNA seeds (default 512bp)
        """
        self.seed_len = seed_len
        self.complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    def augment_sequence(
        self, 
        seq: str, 
        error_rate: Optional[float] = None
    ) -> Tuple[str, str]:
        """
        Three essential augmentations:
        1. Error injection (1-10% - mimics sequencing errors)
        2. Position shift (±50bp - translation continuity)
        3. RC augmentation (50% - strand invariance)
        
        Args:
            seq: DNA sequence string
            error_rate: Error rate (if None, random 1-10%)
        
        Returns:
            (noisy_seq, shifted_seq): Augmented versions
        """
        if error_rate is None:
            error_rate = random.uniform(0.01, 0.10)
        
        # 1. Add realistic sequencing errors
        noisy = self.add_errors(
            seq, 
            error_rate,
            sub_prob=0.6,  # 60% substitutions
            ins_prob=0.2,  # 20% insertions
            del_prob=0.2   # 20% deletions
        )
        
        # 2. Position shift for translation continuity
        # NeuralAligner: ±Lseed/10 exactly (±51bp for 512bp seeds)
        shift = random.randint(-51, 51)
        shifted = self.shift_sequence(seq, shift)
        shifted = self.pad_or_trim(shifted, len(seq))
        
        # 3. Reverse complement (50% chance)
        if random.random() < 0.5:
            shifted = self.reverse_complement(shifted)
        
        return noisy, shifted
    
    def add_errors(
        self, 
        seq: str, 
        error_rate: float,
        sub_prob: float = 0.6,
        ins_prob: float = 0.2,
        del_prob: float = 0.2
    ) -> str:
        """
        Add realistic sequencing errors
        
        Error distribution matches ONT profile:
        - 60% substitutions
        - 20% insertions
        - 20% deletions
        
        Args:
            seq: Input DNA sequence
            error_rate: Overall error rate (0.01-0.10)
            sub_prob: Probability of substitution
            ins_prob: Probability of insertion
            del_prob: Probability of deletion
        
        Returns:
            noisy_seq: Sequence with errors
        """
        assert abs(sub_prob + ins_prob + del_prob - 1.0) < 1e-6, \
            "Error probabilities must sum to 1"
        
        noisy = list(seq)
        i = 0
        
        while i < len(noisy):
            if random.random() < error_rate:
                # Choose error type
                error_type = random.choices(
                    ['sub', 'ins', 'del'],
                    weights=[sub_prob, ins_prob, del_prob]
                )[0]
                
                if error_type == 'sub':
                    # Substitution: replace with different base
                    original = noisy[i]
                    choices = [b for b in ['A', 'C', 'G', 'T'] if b != original]
                    noisy[i] = random.choice(choices)
                    i += 1
                
                elif error_type == 'ins':
                    # Insertion: add random base
                    noisy.insert(i, random.choice(['A', 'C', 'G', 'T']))
                    i += 1  # Skip the inserted base
                    i += 1  # Move to next original position
                
                else:  # deletion
                    # Deletion: remove base
                    if len(noisy) > 1:  # Don't delete if only one base left
                        noisy.pop(i)
                    else:
                        i += 1
            else:
                i += 1
        
        return ''.join(noisy)
    
    def shift_sequence(self, seq: str, shift: int) -> str:
        """
        Shift sequence position (for translation continuity)
        
        Args:
            seq: Input sequence
            shift: Shift amount (negative=left, positive=right)
        
        Returns:
            shifted: Shifted sequence
        """
        if shift > 0:
            # Shift right: remove from start
            return seq[shift:]
        elif shift < 0:
            # Shift left: remove from end
            return seq[:shift]
        else:
            return seq
    
    def pad_or_trim(self, seq: str, target_len: int) -> str:
        """
        Pad or trim sequence to target length
        
        Args:
            seq: Input sequence
            target_len: Target length
        
        Returns:
            adjusted: Sequence of target length
        """
        if len(seq) < target_len:
            # Pad with N's
            padding = 'N' * (target_len - len(seq))
            return seq + padding
        elif len(seq) > target_len:
            # Trim
            return seq[:target_len]
        else:
            return seq
    
    def reverse_complement(self, seq: str) -> str:
        """
        Generate reverse complement
        
        Args:
            seq: Input DNA sequence
        
        Returns:
            rc: Reverse complement
        """
        return ''.join(self.complement.get(base, 'N') for base in reversed(seq))
    
    def create_training_pairs(
        self, 
        genome: str, 
        batch_size: int = 8192
    ) -> Tuple[List[str], List[str]]:
        """
        Generate training pairs for InfoNCE contrastive learning
        
        NO EXPLICIT NEGATIVES! InfoNCE uses other samples in batch as negatives.
        
        Args:
            genome: Reference genome sequence
            batch_size: Number of pairs (large batch crucial - default 8192)
        
        Returns:
            (anchors, positives): Training pairs only
        """
        anchors = []
        positives = []
        
        attempts = 0
        max_attempts = batch_size * 2  # Safety limit
        
        while len(anchors) < batch_size and attempts < max_attempts:
            attempts += 1
            
            # Sample random position for anchor
            pos = random.randint(0, len(genome) - self.seed_len)
            anchor = genome[pos:pos + self.seed_len]
            
            # Skip if contains too many N's
            if anchor.count('N') > self.seed_len * 0.1:
                continue
            
            # Positive: augmented version of same sequence
            error_rate = random.uniform(0.01, 0.10)
            positive = self.add_errors(anchor, error_rate)
            
            # Apply position shift (±51bp for translation continuity)
            shift = random.randint(-51, 51)
            if shift > 0:
                positive = positive[shift:] + 'N' * shift
            elif shift < 0:
                positive = 'N' * abs(shift) + positive[:shift]
            
            # Ensure length is correct
            positive = self.pad_or_trim(positive, self.seed_len)
            
            # RC augmentation (50% for each independently)
            if random.random() < 0.5:
                anchor = self.reverse_complement(anchor)
            if random.random() < 0.5:
                positive = self.reverse_complement(positive)
            
            anchors.append(anchor)
            positives.append(positive)
        
        return anchors, positives
    
    def validate_augmentation(self, seq: str, n_samples: int = 10) -> Dict:
        """
        Validate augmentation quality (spot checks)
        
        Args:
            seq: Test sequence
            n_samples: Number of augmented samples to generate
        
        Returns:
            stats: Statistics about augmentation quality
        """
        stats = {
            'original_len': len(seq),
            'error_rates': [],
            'rc_count': 0,
            'shift_amounts': []
        }
        
        for _ in range(n_samples):
            noisy, shifted = self.augment_sequence(seq)
            
            # Measure error rate
            errors = sum(1 for a, b in zip(seq[:len(noisy)], noisy) if a != b)
            error_rate = errors / len(noisy)
            stats['error_rates'].append(error_rate)
            
            # Check if RC
            if shifted == self.reverse_complement(seq[:len(shifted)]):
                stats['rc_count'] += 1
            
            # Measure shift (approximate)
            shift = len(seq) - len(shifted)
            stats['shift_amounts'].append(shift)
        
        # Calculate averages
        stats['avg_error_rate'] = np.mean(stats['error_rates'])
        stats['rc_fraction'] = stats['rc_count'] / n_samples
        stats['avg_shift'] = np.mean(np.abs(stats['shift_amounts']))
        
        return stats


def load_hard_negative_regions(genome_path: str) -> List[Tuple[int, int]]:
    """
    Load hard negative regions (repetitive, low-complexity)
    
    This should load the 983K regions identified in V3.
    For now, returns placeholder - implement actual loading later.
    
    Args:
        genome_path: Path to reference genome
    
    Returns:
        regions: List of (start, end) tuples
    """
    # TODO: Load actual hard negative regions from V3
    # For now, return empty list (will add random negatives only)
    return []


if __name__ == "__main__":
    # Test augmentation pipeline
    print("Testing GenoCache V4 Augmentation Pipeline\n")
    
    # Test sequence
    test_seq = "ACGTACGTACGTACGT" * 32  # 512bp
    
    augmenter = DataAugmentation(seed_len=512)
    
    # Test single augmentation
    print("1. Testing single augmentation:")
    noisy, shifted = augmenter.augment_sequence(test_seq, error_rate=0.05)
    print(f"   Original length: {len(test_seq)}")
    print(f"   Noisy length: {len(noisy)}")
    print(f"   Shifted length: {len(shifted)}")
    print(f"   Original: {test_seq[:50]}...")
    print(f"   Noisy:    {noisy[:50]}...")
    print(f"   Shifted:  {shifted[:50]}...")
    
    # Validate augmentation
    print("\n2. Validating augmentation (10 samples):")
    stats = augmenter.validate_augmentation(test_seq, n_samples=10)
    print(f"   Average error rate: {stats['avg_error_rate']:.3f}")
    print(f"   RC fraction: {stats['rc_fraction']:.2f}")
    print(f"   Average shift: {stats['avg_shift']:.1f} bp")
    
    # Test pair generation (NO NEGATIVES!)
    print("\n3. Testing pair generation (InfoNCE - no explicit negatives):")
    genome = "ACGTACGTACGTACGT" * 10000  # Mock genome
    anchors, positives = augmenter.create_training_pairs(
        genome, 
        batch_size=10
    )
    print(f"   Generated {len(anchors)} pairs (negatives implicit in batch)")
    print(f"   Anchor sample: {anchors[0][:50]}...")
    print(f"   Positive sample: {positives[0][:50]}...")
    
    print("\n✅ Augmentation pipeline test complete!")
