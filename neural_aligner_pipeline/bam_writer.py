#!/usr/bin/env python3
"""
Phase 4: BAM/SAM output writer
Generates proper alignment files compatible with downstream tools
"""

import pysam
from datetime import datetime

class BAMWriter:
    """Write alignments to BAM/SAM format"""
    
    def __init__(self, output_file, reference_file="GRCh38.fa", mode='wb'):
        """
        Initialize BAM writer
        
        Args:
            output_file: Output BAM/SAM file path
            reference_file: Reference genome FASTA
            mode: 'wb' for BAM, 'w' for SAM
        """
        self.output_file = output_file
        self.reference_file = reference_file
        
        # Create header
        header = self._create_header()
        
        # Open output file
        self.bam_file = pysam.AlignmentFile(output_file, mode, header=header)
        
        print(f"BAM writer initialized: {output_file}")
    
    def _create_header(self):
        """Create SAM/BAM header"""
        from Bio import SeqIO
        
        header = {
            'HD': {'VN': '1.6', 'SO': 'unsorted'},
            'SQ': [],
            'PG': [{
                'ID': 'genocache',
                'PN': 'genocache',
                'VN': '0.1.0',
                'CL': 'genocache align',
                'DS': 'Neural alignment pipeline'
            }]
        }
        
        # Add reference sequences
        print("Building BAM header from reference...")
        for record in SeqIO.parse(self.reference_file, "fasta"):
            header['SQ'].append({
                'SN': record.id,
                'LN': len(record.seq)
            })
        
        print(f"  ✓ Header created with {len(header['SQ'])} sequences")
        return header
    
    def write_alignment(self, read_name, read_seq, alignment, quality_string=None):
        """
        Write alignment to BAM file
        
        Args:
            read_name: Read ID
            read_seq: Read sequence
            alignment: Alignment result dictionary
            quality_string: Quality scores (optional)
        """
        # Create aligned segment with header
        a = pysam.AlignedSegment(self.bam_file.header)
        
        # Basic info
        a.query_name = read_name
        a.query_sequence = read_seq
        a.reference_name = alignment.get('chromosome', 'NC_000022.11')  # Default chr22
        a.reference_start = alignment['position']
        a.mapping_quality = alignment.get('mapq', 60)
        
        # CIGAR string
        if 'cigar' in alignment and alignment['cigar']:
            a.cigarstring = alignment['cigar']
        else:
            # Default: all matches
            a.cigarstring = f"{len(read_seq)}M"
        
        # Flag
        if alignment.get('strand', '+') == '-':
            a.flag = 16  # Reverse strand
        else:
            a.flag = 0   # Forward strand
        
        # Quality scores
        if quality_string:
            a.query_qualities = pysam.qualitystring_to_array(quality_string)
        else:
            # Default quality
            a.query_qualities = pysam.qualitystring_to_array('I' * len(read_seq))
        
        # Optional tags
        if 'score' in alignment:
            a.set_tag('AS', alignment['score'], value_type='i')  # Alignment score
        
        if 'identity' in alignment:
            a.set_tag('ID', alignment['identity'], value_type='f')  # Identity
        
        if 'matches' in alignment:
            a.set_tag('NM', alignment.get('mismatches', 0), value_type='i')  # Edit distance
        
        # Write to file
        self.bam_file.write(a)
    
    def close(self):
        """Close BAM file"""
        self.bam_file.close()
        print(f"BAM file written: {self.output_file}")


class SAMWriter(BAMWriter):
    """Write alignments to SAM format (text)"""
    
    def __init__(self, output_file, reference_file="GRCh38.fa"):
        """Initialize SAM writer"""
        super().__init__(output_file, reference_file, mode='w')


def main():
    """Test BAM writer"""
    print("="*70)
    print("BAM Writer - Phase 4")
    print("="*70)
    
    # Create test BAM file
    output_bam = "test_output.bam"
    
    print("\nCreating BAM writer...")
    writer = BAMWriter(output_bam, reference_file="GRCh38.fa")
    
    # Write test alignment
    print("\nWriting test alignments...")
    
    test_alignments = [
        {
            'read_name': 'read_001',
            'read_seq': 'ACGTACGTACGT' * 10,
            'alignment': {
                'chromosome': 'NC_000022.11',
                'position': 20000000,
                'mapq': 60,
                'cigar': '120M',
                'score': 240,
                'identity': 0.98,
                'matches': 118,
                'mismatches': 2,
                'strand': '+'
            }
        },
        {
            'read_name': 'read_002',
            'read_seq': 'TGCATGCATGCA' * 8,
            'alignment': {
                'chromosome': 'NC_000022.11',
                'position': 30000000,
                'mapq': 40,
                'cigar': '96M',
                'score': 180,
                'identity': 0.95,
                'matches': 91,
                'mismatches': 5,
                'strand': '-'
            }
        }
    ]
    
    for test in test_alignments:
        writer.write_alignment(
            test['read_name'],
            test['read_seq'],
            test['alignment']
        )
        print(f"  ✓ Written: {test['read_name']}")
    
    writer.close()
    
    # Verify BAM file
    print("\nVerifying BAM file...")
    bam = pysam.AlignmentFile(output_bam, 'rb')
    
    count = 0
    for read in bam:
        count += 1
        print(f"  Read {count}: {read.query_name} @ {read.reference_name}:{read.reference_start}")
    
    bam.close()
    
    print(f"\n✓ BAM file verified: {count} alignments")
    
    print("\n" + "="*70)
    print("✓ BAM writer test complete")
    print("="*70)


if __name__ == "__main__":
    main()
