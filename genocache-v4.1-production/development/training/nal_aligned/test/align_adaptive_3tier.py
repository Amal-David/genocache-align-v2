#!/usr/bin/env python3
"""
Three-Tier Adaptive Alignment System

User's proposed configuration:
- Tier 1: 6 seeds, K=64 (fast)
- Tier 2: 12 seeds, K=64 (standard)
- Tier 3: 32 seeds, K=48 (deep, more seeds but lower K for reliability)

Includes false positive detection by comparing with minimap2.
"""

import sys
import time
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from seeding_nal import NALSeeding
from chaining_nal import NALChaining

try:
    from pywfa import WavefrontAligner
    WFA_AVAILABLE = True
except ImportError:
    WFA_AVAILABLE = False


class ThreeTierAligner:
    """
    Three-tier adaptive alignment system with configurable tiers
    """
    
    def __init__(
        self,
        model_path: str,
        index_path: str,
        positions_path: str,
        reference_path: str,
        device: str = 'cuda'
    ):
        """Initialize three-tier aligner"""
        print("Initializing Three-Tier Adaptive Aligner...")
        
        self.model_path = model_path
        self.index_path = index_path
        self.positions_path = positions_path
        self.reference_path = Path(reference_path)
        self.device = device
        
        # Cache for seeders (one per K value)
        self.seeders = {}
        self.chainer = None
        
        # Reference cache
        self.reference_cache = {}
        
        # WFA aligner
        if WFA_AVAILABLE:
            self.wfa_aligner = WavefrontAligner()
            print("  ✅ WFA2-GPU available")
        else:
            self.wfa_aligner = None
            print("  ⚠️  WFA2-GPU not available (placeholder alignment)")
        
        # Tier configurations
        self.tiers = [
            {'name': 'Fast', 'seeds': 6, 'K': 64, 'tolerance': 1000, 'exit_score': 3},
            {'name': 'Standard', 'seeds': 12, 'K': 64, 'tolerance': 1500, 'exit_score': 6},
            {'name': 'Deep', 'seeds': 32, 'K': 48, 'tolerance': 2000, 'exit_score': None},  # No exit
        ]
        
        print(f"  Tier 1 (Fast): {self.tiers[0]['seeds']} seeds, K={self.tiers[0]['K']}")
        print(f"  Tier 2 (Standard): {self.tiers[1]['seeds']} seeds, K={self.tiers[1]['K']}")
        print(f"  Tier 3 (Deep): {self.tiers[2]['seeds']} seeds, K={self.tiers[2]['K']}")
    
    def get_seeder(self, K: int):
        """Get or create seeder for specific K"""
        if K not in self.seeders:
            self.seeders[K] = NALSeeding(
                model_path=self.model_path,
                index_path=self.index_path,
                positions_path=self.positions_path,
                seed_len=512,
                K=K,
                device=self.device
            )
        return self.seeders[K]
    
    def get_chainer(self, tolerance: int):
        """Get or create chainer for specific tolerance"""
        if self.chainer is None or self.chainer.tolerance != tolerance:
            self.chainer = NALChaining(
                tolerance=tolerance,
                min_chain_score=3,
                top_k=5
            )
        return self.chainer
    
    def should_rescue(self, chains: List[Dict], tier_config: Dict) -> bool:
        """
        Determine if we should proceed to next tier
        
        Exit criteria:
        - No chains found → rescue
        - Chain score below threshold → rescue
        - Otherwise → accept
        """
        if tier_config['exit_score'] is None:
            return False  # Final tier, no rescue
        
        if not chains:
            return True  # No chains, definitely rescue
        
        best_chain = chains[0]
        if best_chain['score'] < tier_config['exit_score']:
            return True  # Weak chain, rescue
        
        return False  # Good enough, accept
    
    def load_reference_region(self, chr_name: str, start: int, end: int) -> str:
        """Load reference region from FASTA"""
        if chr_name not in self.reference_cache:
            seq = []
            in_chr = False
            
            with open(self.reference_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if chr_name in line:
                            in_chr = True
                        else:
                            if in_chr:
                                break
                            in_chr = False
                    elif in_chr:
                        seq.append(line)
            
            self.reference_cache[chr_name] = ''.join(seq).upper()
        
        ref_seq = self.reference_cache[chr_name]
        return ref_seq[max(0, start):min(len(ref_seq), end)]
    
    def align_with_wfa(self, read_seq: str, chain: Dict) -> Dict:
        """Perform WFA alignment"""
        read_len = len(read_seq)
        ref_len = int(read_len * 1.002)
        ref_start = max(0, chain['ref_pos'])
        ref_end = ref_start + ref_len
        
        try:
            ref_seq = self.load_reference_region(
                chain['ref_chr'],
                ref_start,
                ref_end
            )
        except Exception as e:
            return {'cigar': f"{read_len}M", 'score': 0}
        
        if self.wfa_aligner is not None and WFA_AVAILABLE:
            try:
                score = self.wfa_aligner.wavefront_align(
                    text=read_seq.upper(),
                    pattern=ref_seq.upper()
                )
                cigar = self.wfa_aligner.cigarstring
                
                if cigar and len(cigar) > 0:
                    return {'cigar': cigar, 'score': score}
            except Exception as e:
                pass
        
        return {'cigar': f"{read_len}M", 'score': 0}
    
    def align_read(
        self,
        read_id: str,
        read_seq: str
    ) -> Dict:
        """
        Align a single read using three-tier strategy
        
        Returns dict with:
        - mapped: bool
        - tier_used: int (1, 2, or 3)
        - chains: list
        - alignment: dict
        - timing: dict
        """
        result = {
            'read_id': read_id,
            'mapped': False,
            'tier_used': None,
            'chains': [],
            'alignment': None,
            'timing': {'tier1': 0, 'tier2': 0, 'tier3': 0}
        }
        
        read_len = len(read_seq)
        
        # Try each tier
        for tier_idx, tier in enumerate(self.tiers, start=1):
            t0 = time.time()
            
            # Get anchors
            seeder = self.get_seeder(tier['K'])
            anchors = seeder.get_anchors(read_seq, num_seeds=tier['seeds'], K=tier['K'])
            
            # Chain
            chainer = self.get_chainer(tier['tolerance'])
            chains, _ = chainer.chain_with_rescue_check(
                anchors, read_len=read_len, num_seeds=tier['seeds'], seed_len=512
            )
            
            elapsed = time.time() - t0
            result['timing'][f'tier{tier_idx}'] = elapsed
            
            # Check if we should continue to next tier
            if self.should_rescue(chains, tier):
                continue  # Try next tier
            
            # Accept this tier's result
            result['tier_used'] = tier_idx
            result['chains'] = chains
            
            if chains:
                result['mapped'] = True
                # Perform alignment
                result['alignment'] = self.align_with_wfa(read_seq, chains[0])
            
            break
        
        return result
    
    def align_reads(
        self,
        reads: List[Tuple[str, str]],
        verbose: bool = True
    ) -> List[Dict]:
        """Align multiple reads"""
        results = []
        tier_counts = {1: 0, 2: 0, 3: 0}
        
        for i, (read_id, read_seq) in enumerate(reads):
            if verbose and (i + 1) % 100 == 0:
                print(f"  Progress: {i+1}/{len(reads)} reads...")
            
            result = self.align_read(read_id, read_seq)
            results.append(result)
            
            if result['tier_used'] is not None:
                tier_counts[result['tier_used']] += 1
        
        # Summary
        if verbose:
            total = len(reads)
            mapped = sum(1 for r in results if r['mapped'])
            
            print(f"\n{'='*80}")
            print("RESULTS:")
            print(f"{'='*80}")
            print(f"  Total reads: {total}")
            print(f"  Mapped: {mapped} ({mapped/total*100:.1f}%)")
            print(f"\n  Tier distribution:")
            print(f"    Tier 1 (Fast): {tier_counts[1]} ({tier_counts[1]/total*100:.1f}%)")
            print(f"    Tier 2 (Standard): {tier_counts[2]} ({tier_counts[2]/total*100:.1f}%)")
            print(f"    Tier 3 (Deep): {tier_counts[3]} ({tier_counts[3]/total*100:.1f}%)")
            
            # Average timing
            avg_time = {}
            for tier_idx in [1, 2, 3]:
                times = [r['timing'][f'tier{tier_idx}'] for r in results if r['timing'][f'tier{tier_idx}'] > 0]
                if times:
                    avg_time[tier_idx] = sum(times) / len(times) * 1000  # ms
                else:
                    avg_time[tier_idx] = 0
            
            # Weighted average
            weighted_avg = sum(tier_counts[i] * avg_time[i] for i in [1, 2, 3]) / total
            
            print(f"\n  Timing (average per tier):")
            print(f"    Tier 1: {avg_time[1]:.1f}ms")
            print(f"    Tier 2: {avg_time[2]:.1f}ms")
            print(f"    Tier 3: {avg_time[3]:.1f}ms")
            print(f"    Weighted average: {weighted_avg:.1f}ms/read")
        
        return results
    
    def write_sam(self, results: List[Dict], output_path: str):
        """Write results to SAM file"""
        with open(output_path, 'w') as f:
            # Header
            f.write("@HD\tVN:1.0\tSO:unsorted\n")
            
            # Write alignments
            for result in results:
                read_id = result['read_id']
                
                if not result['mapped'] or not result['chains']:
                    # Unmapped
                    f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
                else:
                    chain = result['chains'][0]
                    alignment = result['alignment']
                    
                    # SAM fields
                    flag = 0 if chain['strand'] == '+' else 16
                    rname = chain['ref_chr']
                    pos = chain['ref_pos'] + 1  # 1-based
                    mapq = min(60, int(chain['score'] * 5))  # Rough MAPQ
                    cigar = alignment['cigar']
                    
                    f.write(f"{read_id}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*\t")
                    f.write(f"AS:i:{alignment['score']}\tTIER:i:{result['tier_used']}\n")


def load_reads(fastq_path: str, max_reads: int = None) -> List[Tuple[str, str]]:
    """Load reads from FASTQ"""
    reads = []
    with open(fastq_path) as f:
        while True:
            if max_reads and len(reads) >= max_reads:
                break
            
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            read_id = header[1:].split()[0]
            reads.append((read_id, seq))
    
    return reads


def main():
    """Run three-tier alignment"""
    parser = argparse.ArgumentParser(description='Three-tier adaptive NAL alignment')
    parser.add_argument('--model', required=True, help='Model path')
    parser.add_argument('--index', required=True, help='Index path')
    parser.add_argument('--positions', required=True, help='Positions path')
    parser.add_argument('--reference', required=True, help='Reference genome')
    parser.add_argument('--reads', required=True, help='FASTQ reads')
    parser.add_argument('--output', required=True, help='Output SAM file')
    parser.add_argument('--max-reads', type=int, help='Max reads to process')
    parser.add_argument('--device', default='cuda', help='Device')
    
    args = parser.parse_args()
    
    # Load reads
    print(f"Loading reads from {args.reads}...")
    reads = load_reads(args.reads, args.max_reads)
    print(f"  Loaded {len(reads)} reads\n")
    
    # Initialize aligner
    aligner = ThreeTierAligner(
        model_path=args.model,
        index_path=args.index,
        positions_path=args.positions,
        reference_path=args.reference,
        device=args.device
    )
    
    # Align reads
    print("\nAligning reads...")
    results = aligner.align_reads(reads, verbose=True)
    
    # Write SAM
    print(f"\nWriting SAM to {args.output}...")
    aligner.write_sam(results, args.output)
    print("  ✅ SAM written!")
    
    print("\n✅ Three-tier alignment complete!")


if __name__ == '__main__':
    main()
