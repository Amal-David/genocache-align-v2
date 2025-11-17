"""
GenoCache V4 - Complete Pipeline
Adaptive Seeding + WFA Alignment for full SAM/BAM output
"""

import torch
import pickle
import faiss
from Bio import SeqIO
from tqdm import tqdm
import time
from typing import List, Dict, Optional
import argparse

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

from models.encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder
from wfa_alignment import WFAAligner

def load_model(checkpoint_path, device='cuda'):
    """Load trained model from checkpoint"""
    model = GenoCacheEncoder(emb_dim=128, seed_len=512)  # Model was trained with 128D
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(device)
    model.eval()
    return model

class CompletePipeline:
    """
    Complete GenoCache pipeline:
    1. Load model and index
    2. Adaptive seeding (5-16 seeds per read)
    3. Seed chaining
    4. WFA alignment (precise base-level)
    5. SAM/BAM output
    """
    
    def __init__(self,
                 model_path: str,
                 index_path: str,
                 metadata_path: str,
                 genome_path: str,
                 device: str = 'cuda'):
        """
        Args:
            model_path: Path to trained model checkpoint
            index_path: Path to FAISS index
            metadata_path: Path to index metadata
            genome_path: Path to reference genome
            device: 'cuda' or 'cpu'
        """
        print("═" * 80)
        print("GenoCache V4 - Complete Pipeline")
        print("═" * 80)
        print()
        
        # Load model
        print("Loading model...")
        self.model = load_model(model_path, device=device)
        num_params = sum(p.numel() for p in self.model.parameters())
        print(f"✅ Model loaded: {num_params:,} parameters")
        print()
        
        # Load index
        print("Loading FAISS index...")
        self.index = faiss.read_index(index_path)
        # Keep index on CPU for now (GPU FAISS requires faiss-gpu package)
        print(f"✅ Index loaded (CPU mode)")
        print()
        
        # Load metadata
        print("Loading metadata...")
        with open(metadata_path, 'rb') as f:
            self.metadata = pickle.load(f)
        print(f"✅ Metadata loaded:")
        print(f"   Embeddings: {len(self.metadata['positions']):,}")
        print(f"   Chromosomes: {len(set(self.metadata['chr_names']))}")
        print()
        
        # Initialize adaptive seeder
        print("Initializing adaptive seeder...")
        self.seeder = AdaptiveSeeder(
            model=self.model,
            index=self.index,
            metadata=self.metadata,
            min_seeds=5,
            max_seeds=16,
            window_size=512,
            top_k=32
        )
        print("✅ Adaptive seeder ready")
        print()
        
        # Initialize WFA aligner
        print("Initializing WFA aligner...")
        self.aligner = WFAAligner(
            genome_path=genome_path,
            mode='edlib'  # Using edlib for now
        )
        print("✅ WFA aligner ready")
        print()
    
    def align_read(self, read_seq: str, read_id: str) -> Optional[Dict]:
        """
        Align single read through complete pipeline
        
        Args:
            read_seq: Read sequence
            read_id: Read identifier
        
        Returns:
            Complete alignment dict or None if unmapped
        """
        # Step 1: Adaptive seeding
        seed_result = self.seeder.align_read(read_seq, read_id=read_id)
        
        if seed_result is None:
            return None  # Unmapped
        
        # Step 2: WFA precise alignment
        wfa_result = self.aligner.align(
            query=read_seq,
            chr_name=seed_result['chr'],
            start=seed_result['start'],
            end=seed_result['end']
        )
        
        if wfa_result is None:
            return None  # Alignment failed
        
        # Combine results
        result = {
            'read_id': read_id,
            'chr': seed_result['chr'],
            'start': seed_result['start'],
            'end': seed_result['end'],
            'num_seeds': seed_result['num_seeds'],
            'seed_score': seed_result['score'],
            'edit_distance': wfa_result['edit_distance'],
            'identity': wfa_result['identity'],
            'cigar': wfa_result.get('cigar', ''),
            'mapq': int(255 * wfa_result['identity']),
            'status': seed_result['status']
        }
        
        return result
    
    def align_reads(self, 
                   reads_path: str,
                   output_sam: str,
                   max_reads: Optional[int] = None) -> Dict:
        """
        Align all reads from FASTA/FASTQ file
        
        Args:
            reads_path: Path to reads file
            output_sam: Output SAM file path
            max_reads: Maximum reads to process (None = all)
        
        Returns:
            Summary statistics dict
        """
        print("═" * 80)
        print("Aligning reads...")
        print("═" * 80)
        print()
        
        # Load reads
        print(f"Loading reads from {reads_path}...")
        reads = []
        for record in SeqIO.parse(reads_path, "fasta"):
            reads.append((record.id, str(record.seq)))
            if max_reads and len(reads) >= max_reads:
                break
        print(f"✅ Loaded {len(reads)} reads")
        print()
        
        # Process reads
        print("Processing reads...")
        results = []
        unmapped = 0
        
        start_time = time.time()
        
        for read_id, read_seq in tqdm(reads, desc="Aligning"):
            result = self.align_read(read_seq, read_id)
            if result is None:
                unmapped += 1
            else:
                results.append(result)
        
        elapsed = time.time() - start_time
        
        print()
        print("✅ Alignment complete!")
        print(f"   Total reads: {len(reads)}")
        print(f"   Mapped: {len(results)} ({100*len(results)/len(reads):.1f}%)")
        print(f"   Unmapped: {unmapped} ({100*unmapped/len(reads):.1f}%)")
        print(f"   Time: {elapsed:.2f}s")
        print(f"   Speed: {len(reads)/elapsed:.1f} reads/sec")
        print()
        
        # Write SAM file
        print(f"Writing SAM file to {output_sam}...")
        with open(output_sam, 'w') as f:
            # Write SAM header
            f.write("@HD\tVN:1.0\tSO:unsorted\n")
            for chr_name in sorted(set(self.metadata['chr_names'])):
                chr_len = max(self.metadata['positions'][i] 
                            for i, c in enumerate(self.metadata['chr_names']) 
                            if c == chr_name)
                f.write(f"@SQ\tSN:{chr_name}\tLN:{chr_len}\n")
            f.write("@PG\tID:geocache\tPN:GenoCache\tVN:4.0\n")
            
            # Write alignments
            for result in results:
                sam_line = self.aligner.format_sam(
                    read_id=result['read_id'],
                    read_seq='',  # Will need to get from original
                    alignment=result
                )
                f.write(sam_line + "\n")
            
            # Write unmapped reads
            for read_id, read_seq in reads:
                if not any(r['read_id'] == read_id for r in results):
                    f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
        
        print("✅ SAM file written")
        print()
        
        # Summary statistics
        if len(results) > 0:
            avg_identity = sum(r['identity'] for r in results) / len(results)
            avg_seeds = sum(r['num_seeds'] for r in results) / len(results)
            avg_edit = sum(r['edit_distance'] for r in results) / len(results)
            
            print("Summary Statistics:")
            print(f"  Average identity: {100*avg_identity:.1f}%")
            print(f"  Average seeds used: {avg_seeds:.1f}")
            print(f"  Average edit distance: {avg_edit:.1f}")
            print()
        
        return {
            'total_reads': len(reads),
            'mapped': len(results),
            'unmapped': unmapped,
            'time': elapsed,
            'reads_per_sec': len(reads) / elapsed
        }

def main():
    parser = argparse.ArgumentParser(description='GenoCache V4 Complete Pipeline')
    parser.add_argument('--model', type=str, required=True,
                       help='Path to model checkpoint')
    parser.add_argument('--index', type=str, required=True,
                       help='Path to FAISS index')
    parser.add_argument('--metadata', type=str, required=True,
                       help='Path to index metadata')
    parser.add_argument('--genome', type=str, required=True,
                       help='Path to reference genome')
    parser.add_argument('--reads', type=str, required=True,
                       help='Path to reads file (FASTA/FASTQ)')
    parser.add_argument('--output', type=str, required=True,
                       help='Output SAM file')
    parser.add_argument('--max-reads', type=int, default=None,
                       help='Maximum reads to process')
    parser.add_argument('--device', type=str, default='cuda',
                       choices=['cuda', 'cpu'],
                       help='Device to use')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = CompletePipeline(
        model_path=args.model,
        index_path=args.index,
        metadata_path=args.metadata,
        genome_path=args.genome,
        device=args.device
    )
    
    # Align reads
    stats = pipeline.align_reads(
        reads_path=args.reads,
        output_sam=args.output,
        max_reads=args.max_reads
    )
    
    print("═" * 80)
    print("Pipeline complete!")
    print("═" * 80)

if __name__ == '__main__':
    main()
