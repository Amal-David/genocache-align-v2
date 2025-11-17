#!/usr/bin/env python3
"""
Phase 4: Mapping quality (MAPQ) calculation
MAPQ = -10 * log10(P_error)
"""

import numpy as np
import math

class QualityScorer:
    """Calculate mapping quality scores for alignments"""
    
    def __init__(self):
        """Initialize quality scorer"""
        pass
    
    def calculate_mapq(self, primary_alignment, alternative_alignments=None):
        """
        Calculate MAPQ (Mapping Quality) score
        
        MAPQ = -10 * log10(P_error)
        where P_error is probability alignment is incorrect
        
        Args:
            primary_alignment: Best alignment result
            alternative_alignments: List of alternative alignments
        
        Returns:
            MAPQ score (0-60)
        """
        # Start with base quality from alignment identity
        identity = primary_alignment.get('identity', 0.95)
        base_mapq = self._identity_to_mapq(identity)
        
        # Adjust based on SW score
        sw_score = primary_alignment.get('score', 0)
        score_bonus = min(10, sw_score / 1000)  # Up to +10 for high scores
        
        mapq = base_mapq + score_bonus
        
        # Check for multi-mapping (alternative alignments)
        if alternative_alignments and len(alternative_alignments) > 0:
            # Get second-best alignment
            second_best = alternative_alignments[0]
            second_score = second_best.get('score', 0)
            primary_score = sw_score
            
            # Score difference indicates uniqueness
            if primary_score > 0:
                score_ratio = second_score / primary_score if primary_score > 0 else 1.0
                
                if score_ratio > 0.95:
                    # Very similar scores = ambiguous
                    mapq = min(mapq, 3)
                elif score_ratio > 0.90:
                    # Close scores = low confidence
                    mapq = min(mapq, 10)
                elif score_ratio > 0.80:
                    # Moderate difference
                    mapq = min(mapq, 20)
                # else: keep high MAPQ for unique alignment
        
        # Cap at 60 (standard maximum)
        mapq = int(min(60, max(0, mapq)))
        
        return mapq
    
    def _identity_to_mapq(self, identity):
        """
        Convert alignment identity to base MAPQ
        
        Args:
            identity: Alignment identity (0-1)
        
        Returns:
            Base MAPQ score
        """
        if identity >= 0.99:
            return 50
        elif identity >= 0.98:
            return 45
        elif identity >= 0.97:
            return 40
        elif identity >= 0.95:
            return 35
        elif identity >= 0.93:
            return 30
        elif identity >= 0.90:
            return 25
        elif identity >= 0.85:
            return 20
        elif identity >= 0.80:
            return 15
        elif identity >= 0.70:
            return 10
        else:
            return 5
    
    def calculate_alignment_quality(self, alignment):
        """
        Calculate comprehensive quality metrics
        
        Args:
            alignment: Alignment result dictionary
        
        Returns:
            Quality metrics dictionary
        """
        identity = alignment.get('identity', 0)
        matches = alignment.get('matches', 0)
        mismatches = alignment.get('mismatches', 0)
        length = alignment.get('length', 0)
        score = alignment.get('score', 0)
        
        # Calculate various quality metrics
        quality = {
            'identity': identity,
            'matches': matches,
            'mismatches': mismatches,
            'length': length,
            'score': score,
            'match_rate': matches / length if length > 0 else 0,
            'error_rate': mismatches / length if length > 0 else 1.0,
            'quality_category': self._categorize_quality(identity)
        }
        
        return quality
    
    def _categorize_quality(self, identity):
        """Categorize alignment quality"""
        if identity >= 0.99:
            return 'excellent'
        elif identity >= 0.95:
            return 'good'
        elif identity >= 0.90:
            return 'acceptable'
        elif identity >= 0.80:
            return 'poor'
        else:
            return 'very_poor'


def main():
    """Test quality scorer"""
    print("="*70)
    print("Quality Scorer - Phase 4")
    print("="*70)
    
    scorer = QualityScorer()
    
    # Test case 1: High quality unique alignment
    print("\nTest 1: High quality unique alignment")
    alignment1 = {
        'identity': 0.98,
        'score': 5000,
        'matches': 980,
        'mismatches': 20,
        'length': 1000
    }
    mapq1 = scorer.calculate_mapq(alignment1)
    quality1 = scorer.calculate_alignment_quality(alignment1)
    
    print(f"  Identity: {alignment1['identity']*100:.1f}%")
    print(f"  Score: {alignment1['score']}")
    print(f"  MAPQ: {mapq1}")
    print(f"  Quality: {quality1['quality_category']}")
    
    # Test case 2: Ambiguous alignment (multi-mapping)
    print("\nTest 2: Ambiguous alignment (multi-mapping)")
    alignment2 = {
        'identity': 0.95,
        'score': 4000,
        'matches': 950,
        'mismatches': 50,
        'length': 1000
    }
    alternative2 = [
        {'identity': 0.94, 'score': 3900}  # Very close
    ]
    mapq2 = scorer.calculate_mapq(alignment2, alternative2)
    
    print(f"  Primary identity: {alignment2['identity']*100:.1f}%")
    print(f"  Primary score: {alignment2['score']}")
    print(f"  Alternative score: {alternative2[0]['score']}")
    print(f"  Score ratio: {alternative2[0]['score']/alignment2['score']:.3f}")
    print(f"  MAPQ: {mapq2} (low due to ambiguity)")
    
    # Test case 3: Moderate quality
    print("\nTest 3: Moderate quality alignment")
    alignment3 = {
        'identity': 0.85,
        'score': 2000,
        'matches': 850,
        'mismatches': 150,
        'length': 1000
    }
    mapq3 = scorer.calculate_mapq(alignment3)
    quality3 = scorer.calculate_alignment_quality(alignment3)
    
    print(f"  Identity: {alignment3['identity']*100:.1f}%")
    print(f"  Score: {alignment3['score']}")
    print(f"  MAPQ: {mapq3}")
    print(f"  Quality: {quality3['quality_category']}")
    
    print("\n" + "="*70)
    print("✓ Quality scorer test complete")
    print("="*70)


if __name__ == "__main__":
    main()
