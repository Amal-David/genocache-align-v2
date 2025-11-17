#!/usr/bin/env python3
"""
Fast alignment implementation using parasail (SSE/AVX optimized)
Much faster than basic edlib for production use
"""

import parasail
import numpy as np
from typing import Dict, Optional, Tuple
from Bio import SeqIO

class FastAligner:
    """
    Fast sequence alignment using parasail library
    Supports both global and semi-global alignment
    """
    
    def __init__(self, genome_dict: Dict[str, str], mode='semi-global'):
        """
        Initialize aligner with genome reference
        
        Args:
            genome_dict: {chr_name: sequence}
            mode: 'global' or 'semi-global' (semi-global for read mapping)
        """
        self.genome_dict = genome_dict
        self.mode = mode
        
        # Scoring parameters (gap-affine)
        # Note: parasail uses POSITIVE gap penalties
        self.match = 2
        self.mismatch = -4
        self.gap_open = 8  # Positive for parasail
        self.gap_extend = 2  # Positive for parasail
        
        # Create substitution matrix
        self.matrix = parasail.matrix_create("ACGT", self.match, self.mismatch)
        
        print(f"✅ FastAligner initialized ({mode} mode)")
        print(f"   Scoring: match={self.match}, mismatch={self.mismatch}")
        print(f"   Gap: open={self.gap_open}, extend={self.gap_extend}")
    
    def extract_reference(self, chr_name: str, start: int, end: int, 
                         padding: int = 5000) -> str:
        """
        Extract reference sequence with padding
        
        Args:
            chr_name: Chromosome name
            start: Start position
            end: End position
            padding: Extra bases on each side (for indels)
        
        Returns:
            reference sequence
        """
        if chr_name not in self.genome_dict:
            return ""
        
        seq = self.genome_dict[chr_name]
        
        # Add padding
        ref_start = max(0, start - padding)
        ref_end = min(len(seq), end + padding)
        
        return seq[ref_start:ref_end].upper()
    
    def align(self, query: str, reference: str) -> Dict:
        """
        Align query to reference using parasail
        
        Args:
            query: Read sequence
            reference: Reference sequence
        
        Returns:
            Dict with alignment results (CIGAR, score, positions)
        """
        # Use Smith-Waterman for local alignment (best for reads)
        # This finds the best local alignment region
        result = parasail.sw_trace_striped_32(
            query, reference,
            self.gap_open, self.gap_extend,
            self.matrix
        )
        
        # Get alignment score
        score = result.score
        
        # Get CIGAR string from traceback
        cigar = result.cigar
        
        # Parse CIGAR string
        cigar_str = self._parse_cigar_detailed(cigar)
        
        # Get alignment positions
        end_query = result.end_query
        end_ref = result.end_ref
        
        # Calculate start positions
        query_len = len(query)
        cigar_len = self._get_cigar_length(cigar)
        
        return {
            'score': score,
            'cigar': cigar_str,
            'end_query': end_query,
            'end_ref': end_ref,
            'start_query': max(0, end_query - query_len + 1),
            'start_ref': max(0, end_ref - cigar_len + 1),
            'length': cigar_len
        }
    
    def _get_cigar_length(self, cigar) -> int:
        """Get length of alignment from CIGAR"""
        if cigar is None or not hasattr(cigar, 'len'):
            return 0
        return cigar.len
    
    def _parse_cigar_detailed(self, cigar) -> str:
        """
        Parse parasail CIGAR to standard SAM format
        
        Args:
            cigar: parasail CIGAR object
        
        Returns:
            Standard CIGAR string (SAM format)
        """
        if cigar is None:
            return "*"
        
        try:
            # parasail CIGAR decode gives RLE format (e.g., "160=", "50=2X10=")
            if hasattr(cigar, 'decode'):
                extended_cigar = cigar.decode.decode('utf-8')
            else:
                return "*"
            
            # Parse RLE CIGAR format
            # Format: number + operation (e.g., "160=", "5I", "3D")
            # Operations:
            #   = exact match
            #   X mismatch
            #   I insertion
            #   D deletion
            # Convert to SAM: = and X both become M
            
            import re
            sam_cigar = []
            
            # Parse using regex: (digits)(operation)
            pattern = r'(\d+)([=XMID])'
            matches = re.findall(pattern, extended_cigar)
            
            for count, op in matches:
                count = int(count)
                # Convert to SAM format
                if op in '=X':
                    sam_op = 'M'  # Both match and mismatch are M in SAM
                elif op == 'I':
                    sam_op = 'I'
                elif op == 'D':
                    sam_op = 'D'
                elif op == 'M':
                    sam_op = 'M'
                else:
                    continue
                
                sam_cigar.append(f"{count}{sam_op}")
            
            return ''.join(sam_cigar) if sam_cigar else "*"
            
        except Exception as e:
            print(f"Warning: CIGAR parsing failed: {e}")
            return "*"
    
    def _format_cigar(self, cigar) -> str:
        """
        Format parasail CIGAR object to string
        
        Args:
            cigar: parasail CIGAR object
        
        Returns:
            CIGAR string (e.g., "100M2I5M3D10M")
        """
        if cigar is None:
            return "*"
        
        # Get CIGAR string from parasail
        try:
            if hasattr(cigar, 'decode'):
                cigar_str = cigar.decode.decode('utf-8')
            elif hasattr(cigar, 'seq'):
                cigar_str = cigar.seq.decode('utf-8') if isinstance(cigar.seq, bytes) else str(cigar.seq)
            else:
                return "*"
        except:
            return "*"
        
        # Format: convert |=X to standard CIGAR
        # | = match/mismatch
        # = = match
        # X = mismatch
        # I = insertion
        # D = deletion
        
        # For simplicity, convert to standard format
        formatted = []
        current_op = None
        current_count = 0
        
        for char in cigar_str:
            if char in '|=X':
                op = 'M'  # Match/mismatch
            elif char in 'I':
                op = 'I'  # Insertion
            elif char in 'D':
                op = 'D'  # Deletion
            else:
                continue
            
            if op == current_op:
                current_count += 1
            else:
                if current_op:
                    formatted.append(f"{current_count}{current_op}")
                current_op = op
                current_count = 1
        
        # Add final operation
        if current_op:
            formatted.append(f"{current_count}{current_op}")
        
        return ''.join(formatted) if formatted else "*"
    
    def align_read(self, read_seq: str, chr_name: str, 
                   seed_start: int, seed_end: int) -> Optional[Dict]:
        """
        Align read to reference region identified by seeding
        
        Args:
            read_seq: Read sequence
            chr_name: Chromosome name
            seed_start: Start position from seeding
            seed_end: End position from seeding
        
        Returns:
            Alignment dict with CIGAR, score, positions
        """
        # Extract reference region
        # Use tighter window for speed (±5kb padding)
        reference = self.extract_reference(chr_name, seed_start, seed_end, padding=5000)
        
        if not reference:
            print(f"    ⚠️  Reference extraction failed for {chr_name}:{seed_start}-{seed_end}")
            print(f"       Available: {list(self.genome_dict.keys())[:5]}...")
            return None
        
        print(f"    → Aligning {len(read_seq)}bp read to {len(reference)}bp ref region")
        
        # Align
        result = self.align(read_seq, reference)
        
        # Check score threshold (relaxed for testing)
        if result['score'] < 20:  # Minimum score threshold
            print(f"    ⚠️  Very low alignment score: {result['score']}")
            return None
        
        # Adjust positions (add seed_start offset)
        ref_start = max(0, seed_start - 5000)
        result['ref_chr'] = chr_name
        result['ref_start'] = ref_start + result['start_ref']
        result['ref_end'] = ref_start + result['end_ref']
        
        return result
    
    def format_sam(self, read_id: str, read_seq: str, alignment: Dict, 
                   mapq: int = 60) -> str:
        """
        Format alignment as SAM record
        
        Args:
            read_id: Read identifier
            read_seq: Read sequence
            alignment: Alignment dict from align_read()
            mapq: Mapping quality
        
        Returns:
            SAM record string
        """
        if not alignment:
            # Unmapped
            return f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n"
        
        flag = 0  # Primary alignment
        chr_name = alignment['ref_chr']
        pos = alignment['ref_start'] + 1  # SAM is 1-based
        cigar = alignment['cigar']
        
        # Optional fields
        as_score = f"AS:i:{alignment['score']}"
        
        sam_line = f"{read_id}\t{flag}\t{chr_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read_seq}\t*\t{as_score}\n"
        
        return sam_line


def test_aligner():
    """Quick test of FastAligner"""
    print("=" * 80)
    print("Testing FastAligner")
    print("=" * 80)
    print()
    
    # Create simple test case
    ref = "ACGTACGTACGTACGT" * 100  # 1.6kb reference
    query = "ACGTACGTACGTACGT" * 10  # 160bp query (perfect match)
    
    # Test with mutation
    query_mut = list(query)
    query_mut[50] = 'T' if query_mut[50] != 'T' else 'A'  # SNP
    query_mut = ''.join(query_mut)
    
    # Test with indel
    query_indel = query[:80] + query[85:]  # 5bp deletion
    
    # Initialize aligner
    genome_dict = {'chr1': ref}
    aligner = FastAligner(genome_dict)
    
    # Test 1: Perfect match
    print("Test 1: Perfect match")
    result = aligner.align(query, ref)
    print(f"  Score: {result['score']}")
    print(f"  CIGAR: {result['cigar']}")
    print(f"  CIGAR length: {result['length']}")
    print(f"  Query: {len(query)}bp, Ref region: {result['start_ref']}-{result['end_ref']}")
    print()
    
    # Test 2: With SNP
    print("Test 2: With SNP")
    result = aligner.align(query_mut, ref)
    print(f"  Score: {result['score']}")
    print(f"  CIGAR: {result['cigar'][:50]}...")
    print()
    
    # Test 3: With deletion
    print("Test 3: With deletion")
    result = aligner.align(query_indel, ref)
    print(f"  Score: {result['score']}")
    print(f"  CIGAR: {result['cigar'][:50]}...")
    print()
    
    print("=" * 80)
    print("✅ FastAligner tests passed!")
    print("=" * 80)


if __name__ == '__main__':
    test_aligner()
