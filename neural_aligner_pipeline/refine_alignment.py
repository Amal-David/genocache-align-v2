#!/usr/bin/env python3
"""
Phase 3: Smith-Waterman refinement for precise alignment
Expected improvement: 73% → 90-95%
"""

import parasail
from Bio import SeqIO
import numpy as np

class AlignmentRefiner:
    def __init__(self, ref_fa="GRCh38.fa"):
        """
        Initialize alignment refiner with reference genome
        
        Args:
            ref_fa: Path to reference FASTA file
        """
        print("Loading reference genome for SW refinement...")
        self.reference = {}
        self.chr_names = []
        
        for record in SeqIO.parse(ref_fa, "fasta"):
            chr_name = record.id
            self.reference[chr_name] = str(record.seq).upper()
            self.chr_names.append(chr_name)
        
        print(f"  ✓ Loaded {len(self.reference)} chromosomes")
        
        # For chr22-only testing, find chr22
        self.chr22_name = None
        for name in self.chr_names:
            if "NC_000022" in name or name == "chr22":
                self.chr22_name = name
                print(f"  ✓ Chr22 identified as: {name}")
                break
    
    def refine_alignment(self, read_seq, seed_position, chromosome=None, margin=2000):
        """
        Refine alignment using Smith-Waterman
        
        Args:
            read_seq: Query sequence
            seed_position: Approximate position from seeding
            chromosome: Chromosome name (default: chr22)
            margin: Bases to include around seed position
        
        Returns:
            Refined alignment with CIGAR, position, score
        """
        # Default to chr22 for testing
        if chromosome is None:
            chromosome = self.chr22_name
        
        if chromosome not in self.reference:
            return None
        
        ref_seq = self.reference[chromosome]
        read_len = len(read_seq)
        
        # Extract reference region with margin
        ref_start = max(0, seed_position - margin)
        ref_end = min(len(ref_seq), seed_position + read_len + margin)
        ref_region = ref_seq[ref_start:ref_end]
        
        # Smith-Waterman alignment using parasail
        # Use simpler match/mismatch scoring
        result = parasail.sw_stats_scan_16(
            read_seq,
            ref_region,
            open=3,      # Gap open penalty
            extend=1,    # Gap extension penalty
            matrix=parasail.matrix_create("ACGT", 2, -1)  # +2 match, -1 mismatch
        )
        
        # Calculate actual position in reference
        # SW finds local alignment - end_ref tells us where it ends
        # Alignment starts at: end_ref - alignment_length
        alignment_length = result.length if result.length > 0 else read_len
        
        # Position where alignment starts in reference region
        ref_align_start = result.end_ref - alignment_length + 1
        
        # Convert to global reference position
        aligned_pos = ref_start + ref_align_start
        
        # Calculate identity
        if result.matches > 0 and alignment_length > 0:
            identity = result.matches / alignment_length
        else:
            # Estimate from score
            identity = min(1.0, max(0.0, result.score / (read_len * 2)))
        
        return {
            'position': aligned_pos,
            'cigar': f"{alignment_length}M",  # Simple CIGAR for now
            'score': result.score,
            'matches': result.matches,
            'mismatches': alignment_length - result.matches if alignment_length > result.matches else 0,
            'length': alignment_length,
            'identity': identity,
            'end_query': result.end_query,
            'end_ref': result.end_ref
        }
    
    def result_to_cigar(self, result):
        """
        Convert parasail result to CIGAR string
        
        Args:
            result: Parasail alignment result
        
        Returns:
            CIGAR string (e.g., "100M2D50M1I47M")
        """
        if not hasattr(result, 'cigar') or result.cigar is None:
            # Fallback: simple match CIGAR
            return f"{result.len_query}M"
        
        try:
            cigar_list = []
            cigar = result.cigar
            
            # parasail cigar is encoded, decode it
            for i in range(len(cigar.seq)):
                length = cigar.seq[i]
                op_code = cigar.decode[i]
                
                if op_code == '=':
                    op_code = 'M'  # Convert = to M for SAM compatibility
                elif op_code == 'X':
                    op_code = 'M'  # Convert X to M
                
                cigar_list.append(f"{length}{op_code}")
            
            return ''.join(cigar_list)
        except:
            # Fallback on error
            return f"{result.len_query}M"


class RefinedAligner:
    """Complete aligner with multi-seeding + chaining + SW refinement"""
    
    def __init__(self, chained_aligner, refiner):
        """
        Initialize refined aligner
        
        Args:
            chained_aligner: ChainedMultiSeedAligner instance
            refiner: AlignmentRefiner instance
        """
        self.seeder = chained_aligner
        self.refiner = refiner
    
    def align_read(self, read_seq, n_seeds=5):
        """
        Align read using full pipeline: seed → chain → refine
        
        Args:
            read_seq: Read sequence
            n_seeds: Number of seeds for initial seeding
        
        Returns:
            Refined alignment with exact position and CIGAR
        """
        # Step 1: Seed and chain
        seed_result = self.seeder.align_read(read_seq, n_seeds=n_seeds)
        
        if not seed_result:
            return None
        
        # Step 2: Refine with Smith-Waterman
        refined = self.refiner.refine_alignment(
            read_seq,
            seed_position=seed_result['ref_pos'],
            chromosome=None,  # Use default chr22
            margin=2000
        )
        
        if not refined:
            # Fallback to seed result
            return {
                'position': seed_result['ref_pos'],
                'cigar': None,
                'score': seed_result['chain_score'],
                'chain_length': seed_result['chain_length'],
                'identity': seed_result['avg_score'],
                'refined': False
            }
        
        # Combine results
        return {
            'position': refined['position'],
            'cigar': refined['cigar'],
            'score': refined['score'],
            'matches': refined['matches'],
            'mismatches': refined['mismatches'],
            'identity': refined['identity'],
            'chain_length': seed_result['chain_length'],
            'seed_score': seed_result['chain_score'],
            'refined': True
        }


def main():
    """Test Smith-Waterman refinement"""
    print("="*70)
    print("Smith-Waterman Refinement - Phase 3")
    print("="*70)
    
    # Load refiner
    refiner = AlignmentRefiner(ref_fa="GRCh38.fa")
    
    # Test with synthetic read from chr22
    print("\nExtracting test read from chr22...")
    chr22_seq = refiner.reference[refiner.chr22_name]
    
    # Find a region without too many Ns
    test_pos = 20000000  # Different region
    while test_pos < len(chr22_seq) - 10000:
        test_read = chr22_seq[test_pos:test_pos + 1000]
        if test_read.count('N') < 100:  # Less than 10% Ns
            break
        test_pos += 10000
    
    print(f"  True position: {test_pos:,}")
    print(f"  Read length: {len(test_read)}")
    
    # Add some errors
    import random
    read_list = list(test_read)
    bases = ['A', 'C', 'G', 'T']
    for i in range(len(read_list)):
        if random.random() < 0.05:  # 5% error rate
            read_list[i] = random.choice([b for b in bases if b != read_list[i]])
    test_read_with_errors = ''.join(read_list)
    
    # Test refinement with approximate position (off by 50bp)
    seed_pos = test_pos + 50
    print(f"\nRefining from seed position: {seed_pos:,} (50bp off)")
    
    result = refiner.refine_alignment(test_read_with_errors, seed_pos, margin=2000)
    
    if result:
        error = abs(result['position'] - test_pos)
        print(f"\nRefinement results:")
        print(f"  Refined position: {result['position']:,}")
        print(f"  Error: {error} bp")
        print(f"  SW score: {result['score']}")
        print(f"  Matches: {result['matches']}")
        print(f"  Mismatches: {result['mismatches']}")
        print(f"  Identity: {result['identity']*100:.1f}%")
        print(f"  CIGAR: {result['cigar'][:50]}...")
        
        if error < 10:
            print("\n  ✓ Refinement successful (error < 10bp)")
        elif error < 100:
            print("\n  ~ Refinement good (error < 100bp)")
        else:
            print("\n  ✗ Refinement needs tuning")
    else:
        print("\n  ✗ Refinement failed")
    
    print("\n" + "="*70)
    print("✓ SW refinement test complete")
    print("="*70)


if __name__ == "__main__":
    main()
