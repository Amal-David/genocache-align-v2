#!/usr/bin/env python3
"""
Complete Alignment Pipeline - Phase 4
Multi-seeding → Chaining → SW Refinement → Quality Scoring → BAM Output
"""

from multi_seeder import MultiSeedAligner
from seed_chainer import ChainedMultiSeedAligner
from refine_alignment import AlignmentRefiner, RefinedAligner
from quality_scorer import QualityScorer
from bam_writer import BAMWriter
from Bio import SeqIO
import argparse

class CompleteAligner:
    """
    Complete alignment pipeline with all components
    """
    
    def __init__(self, encoder_path, index_path, positions_path, ref_fa):
        """
        Initialize complete aligner
        
        Args:
            encoder_path: Path to trained encoder model
            index_path: Path to FAISS index
            positions_path: Path to position map
            ref_fa: Path to reference genome
        """
        print("Initializing complete alignment pipeline...")
        
        # Stage 1: Multi-seed aligner
        print("  [1/4] Loading multi-seed aligner...")
        base_aligner = MultiSeedAligner(encoder_path, index_path, positions_path)
        
        # Stage 2: Chained aligner
        print("  [2/4] Initializing chaining...")
        chained_aligner = ChainedMultiSeedAligner(base_aligner, max_gap=15000, max_deviation=2000)
        
        # Stage 3: SW refiner
        print("  [3/4] Loading SW refiner...")
        refiner = AlignmentRefiner(ref_fa)
        self.refined_aligner = RefinedAligner(chained_aligner, refiner)
        
        # Stage 4: Quality scorer
        print("  [4/4] Initializing quality scorer...")
        self.quality_scorer = QualityScorer()
        
        print("✓ Complete pipeline ready")
    
    def align_read(self, read_name, read_seq, n_seeds=5):
        """
        Align a single read through complete pipeline
        
        Args:
            read_name: Read identifier
            read_seq: Read sequence
            n_seeds: Number of seeds to use
        
        Returns:
            Complete alignment with quality scores
        """
        # Get refined alignment
        alignment = self.refined_aligner.align_read(read_seq, n_seeds=n_seeds)
        
        if not alignment:
            return None
        
        # Calculate MAPQ
        mapq = self.quality_scorer.calculate_mapq(alignment)
        alignment['mapq'] = mapq
        
        # Add read info
        alignment['read_name'] = read_name
        alignment['read_seq'] = read_seq
        
        return alignment
    
    def align_reads(self, reads_file, output_bam, n_seeds=5, max_reads=None):
        """
        Align reads from FASTA/FASTQ file
        
        Args:
            reads_file: Input FASTA/FASTQ file
            output_bam: Output BAM file
            n_seeds: Number of seeds per read
            max_reads: Maximum reads to process (None = all)
        
        Returns:
            Alignment statistics
        """
        print(f"\nAligning reads: {reads_file} → {output_bam}")
        
        # Open BAM writer
        writer = BAMWriter(output_bam, reference_file="GRCh38.fa")
        
        # Statistics
        stats = {
            'total': 0,
            'aligned': 0,
            'unaligned': 0,
            'mapq_distribution': {0: 0, 10: 0, 20: 0, 30: 0, 40: 0, 50: 0, 60: 0}
        }
        
        # Process reads
        print("\nProcessing reads...")
        
        for i, record in enumerate(SeqIO.parse(reads_file, "fasta")):
            if max_reads and i >= max_reads:
                break
            
            read_name = record.id
            read_seq = str(record.seq).upper()
            
            # Align
            alignment = self.align_read(read_name, read_seq, n_seeds=n_seeds)
            
            if alignment:
                # Write to BAM
                writer.write_alignment(read_name, read_seq, alignment)
                stats['aligned'] += 1
                
                # Update MAPQ distribution
                mapq = alignment['mapq']
                mapq_bin = (mapq // 10) * 10
                stats['mapq_distribution'][mapq_bin] = stats['mapq_distribution'].get(mapq_bin, 0) + 1
            else:
                stats['unaligned'] += 1
            
            stats['total'] += 1
            
            # Progress
            if (i + 1) % 100 == 0:
                pct_aligned = stats['aligned'] / stats['total'] * 100
                print(f"  Processed: {i+1:6d} | Aligned: {stats['aligned']:6d} ({pct_aligned:5.1f}%)")
        
        writer.close()
        
        # Print statistics
        print(f"\n{'='*70}")
        print("Alignment Statistics:")
        print(f"{'='*70}")
        print(f"  Total reads:    {stats['total']:,}")
        print(f"  Aligned:        {stats['aligned']:,} ({stats['aligned']/stats['total']*100:.1f}%)")
        print(f"  Unaligned:      {stats['unaligned']:,} ({stats['unaligned']/stats['total']*100:.1f}%)")
        print(f"\n  MAPQ distribution:")
        for mapq in sorted(stats['mapq_distribution'].keys()):
            count = stats['mapq_distribution'][mapq]
            if count > 0:
                pct = count / stats['aligned'] * 100 if stats['aligned'] > 0 else 0
                print(f"    {mapq:2d}-{mapq+9:2d}: {count:6d} ({pct:5.1f}%)")
        print(f"{'='*70}")
        
        return stats


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description='Complete alignment pipeline')
    parser.add_argument('--reads', required=True, help='Input FASTA/FASTQ file')
    parser.add_argument('--output', required=True, help='Output BAM file')
    parser.add_argument('--encoder', default='nal_encoder_best.pt', help='Encoder model')
    parser.add_argument('--index', default='faiss_index_improved.idx', help='FAISS index')
    parser.add_argument('--positions', default='ref_positions_improved.npy', help='Position map')
    parser.add_argument('--reference', default='GRCh38.fa', help='Reference genome')
    parser.add_argument('--seeds', type=int, default=5, help='Seeds per read')
    parser.add_argument('--max-reads', type=int, help='Maximum reads to process')
    
    args = parser.parse_args()
    
    # Create aligner
    aligner = CompleteAligner(
        encoder_path=args.encoder,
        index_path=args.index,
        positions_path=args.positions,
        ref_fa=args.reference
    )
    
    # Align reads
    stats = aligner.align_reads(
        reads_file=args.reads,
        output_bam=args.output,
        n_seeds=args.seeds,
        max_reads=args.max_reads
    )
    
    print(f"\n✓ Complete! Output: {args.output}")


if __name__ == "__main__":
    main()
