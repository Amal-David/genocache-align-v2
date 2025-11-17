"""
GenoCache V4 - WFA (Wavefront Alignment) Integration
For precise base-level alignment after seeding
"""

import subprocess
import os
from typing import Dict, Optional, Tuple
from Bio import SeqIO
import tempfile

class WFAAligner:
    """
    Wrapper for Wavefront Alignment Algorithm (WFA)
    
    WFA is exact alignment algorithm with O(s²) complexity where s = edit distance
    Much faster than Smith-Waterman O(nm) for high-identity sequences (ONT)
    
    Note: For now using edlib as fallback (similar performance, easier to install)
    TODO: Integrate actual WFA-GPU when available
    """
    
    def __init__(self, 
                 genome_path: str,
                 mode: str = 'edlib',  # 'edlib', 'wfa', or 'wfa-gpu'
                 k: int = 8):
        """
        Args:
            genome_path: Path to reference genome FASTA
            mode: Alignment mode ('edlib' for now, 'wfa-gpu' later)
            k: k-mer size for seeding (if needed)
        """
        self.genome_path = genome_path
        self.mode = mode
        self.k = k
        
        # Load genome into memory
        print(f"Loading genome from {genome_path}...")
        self.genome = {}
        for record in SeqIO.parse(genome_path, "fasta"):
            self.genome[record.id] = str(record.seq).upper()
        print(f"Loaded {len(self.genome)} chromosomes")
        
        # Check if edlib is available
        try:
            import edlib
            self.edlib = edlib
            self.has_edlib = True
        except ImportError:
            print("Warning: edlib not installed. Install with: pip install edlib")
            self.has_edlib = False
    
    def extract_reference(self, 
                         chr_name: str, 
                         start: int, 
                         end: int,
                         buffer: int = 1000) -> Optional[str]:
        """
        Extract reference sequence for alignment
        
        Args:
            chr_name: Chromosome name
            start: Start position
            end: End position
            buffer: Buffer on each side (default 1kb)
        
        Returns:
            Reference sequence or None if invalid
        """
        if chr_name not in self.genome:
            return None
        
        seq = self.genome[chr_name]
        
        # Add buffer
        ref_start = max(0, start - buffer)
        ref_end = min(len(seq), end + buffer)
        
        return seq[ref_start:ref_end]
    
    def align_edlib(self, query: str, reference: str) -> Dict:
        """
        Align using edlib library (fast CPU alignment)
        
        Args:
            query: Query sequence (read)
            reference: Reference sequence
        
        Returns:
            Alignment dict with edit_distance, cigar, identity
        """
        if not self.has_edlib:
            return {
                'edit_distance': -1,
                'cigar': '',
                'identity': 0.0,
                'error': 'edlib not available'
            }
        
        # Run edlib alignment
        result = self.edlib.align(query, reference, task="path")
        
        if result['editDistance'] == -1:
            return {
                'edit_distance': -1,
                'cigar': '',
                'identity': 0.0,
                'error': 'alignment failed'
            }
        
        # Calculate identity
        edit_distance = result['editDistance']
        identity = 1.0 - (edit_distance / len(query))
        
        # Get CIGAR string
        cigar = self.edlib.getNiceAlignment(result, query, reference)
        
        return {
            'edit_distance': edit_distance,
            'cigar': result.get('cigar', ''),
            'identity': identity,
            'alignment': cigar
        }
    
    def align_wfa(self, query: str, reference: str) -> Dict:
        """
        Align using WFA (Wavefront Alignment)
        
        TODO: Implement actual WFA integration
        For now, falls back to edlib
        """
        return self.align_edlib(query, reference)
    
    def align_wfa_gpu(self, queries: list, references: list) -> list:
        """
        Batch align using WFA-GPU
        
        TODO: Implement WFA-GPU batch alignment
        For now, falls back to edlib
        """
        results = []
        for query, reference in zip(queries, references):
            results.append(self.align_edlib(query, reference))
        return results
    
    def align(self, 
             query: str,
             chr_name: str,
             start: int,
             end: int) -> Optional[Dict]:
        """
        Align query to reference region
        
        Args:
            query: Query sequence (read)
            chr_name: Chromosome name
            start: Start position (from seeding)
            end: End position (from seeding)
        
        Returns:
            Alignment result dict or None
        """
        # Extract reference
        reference = self.extract_reference(chr_name, start, end)
        if reference is None:
            return None
        
        # Align based on mode
        if self.mode == 'wfa-gpu':
            result = self.align_wfa_gpu([query], [reference])[0]
        elif self.mode == 'wfa':
            result = self.align_wfa(query, reference)
        else:  # edlib
            result = self.align_edlib(query, reference)
        
        # Add position info
        result['chr'] = chr_name
        result['start'] = start
        result['end'] = end
        
        return result
    
    def format_sam(self, 
                   read_id: str,
                   read_seq: str,
                   alignment: Dict) -> str:
        """
        Format alignment as SAM record
        
        Args:
            read_id: Read identifier
            read_seq: Read sequence
            alignment: Alignment dict
        
        Returns:
            SAM format string
        """
        if alignment is None or alignment['edit_distance'] == -1:
            # Unmapped
            return f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*"
        
        # SAM fields
        qname = read_id
        flag = 0  # Mapped
        rname = alignment['chr']
        pos = alignment['start'] + 1  # 1-based
        mapq = int(255 * alignment['identity'])  # Simple MAPQ
        cigar = alignment.get('cigar', f"{len(read_seq)}M")
        rnext = "*"
        pnext = 0
        tlen = 0
        seq = read_seq
        qual = "*"
        
        # Optional fields
        tags = f"NM:i:{alignment['edit_distance']}\tID:f:{alignment['identity']:.4f}"
        
        sam_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{tags}"
        
        return sam_line

class WFAGPUBatcher:
    """
    Batch processor for WFA-GPU
    
    Processes multiple alignments in parallel on GPU
    TODO: Implement when WFA-GPU is available
    """
    
    def __init__(self, batch_size: int = 256):
        self.batch_size = batch_size
    
    def align_batch(self, queries: list, references: list) -> list:
        """
        Batch align on GPU
        
        Args:
            queries: List of query sequences
            references: List of reference sequences
        
        Returns:
            List of alignment results
        """
        # TODO: Implement WFA-GPU batch processing
        # For now, use edlib in batches
        results = []
        for query, reference in zip(queries, references):
            # Placeholder
            results.append({
                'edit_distance': -1,
                'cigar': '',
                'identity': 0.0
            })
        return results

if __name__ == '__main__':
    print("WFAAligner module ready")
    print("Modes: 'edlib' (current), 'wfa', 'wfa-gpu' (future)")
    print()
    print("Usage:")
    print("  aligner = WFAAligner(genome_path='GRCh38.fa')")
    print("  result = aligner.align(query, chr_name, start, end)")
    print("  sam = aligner.format_sam(read_id, read_seq, result)")
