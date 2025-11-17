#!/usr/bin/env python3
"""
GenoCache Alignment - minimap2-Compatible Output

Enhanced version with full SAM tags and secondary alignments
"""

import sys
import argparse
import torch
import faiss
import pickle
from pathlib import Path
from Bio import SeqIO
import time

# Add core module to path
sys.path.insert(0, str(Path(__file__).parent / "genocache_core"))

from encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner
from sam_output import write_sam_with_secondaries, write_sam_header

def load_model(model_path, device):
    """Load GenoCache encoder model"""
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    checkpoint = torch.load(model_path, map_location=device)
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model.load_state_dict(checkpoint)
    
    model.eval()
    return model


def main():
    parser = argparse.ArgumentParser(description='GenoCache Alignment - minimap2 Compatible')
    parser.add_argument('--reads', required=True, help='Input FASTQ file')
    parser.add_argument('--output', required=True, help='Output SAM file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('--model', default='models/genocache_model.pt', help='Model path')
    parser.add_argument('--index', default='indexes/genocache_v4_production.index', help='Index path')
    parser.add_argument('--metadata', default='indexes/genocache_v4_production.metadata.pkl', help='Metadata path')
    parser.add_argument('--secondary', action='store_true', help='Output secondary alignments')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'], help='Device')
    
    args = parser.parse_args()
    
    print("="*80)
    print("GenoCache Alignment - minimap2 Compatible")
    print("="*80)
    print()
    print(f"Reads:      {args.reads}")
    print(f"Reference:  {args.reference}")
    print(f"Output:     {args.output}")
    print(f"Secondary:  {args.secondary}")
    print()
    
    # Load model
    print("Loading model...")
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    model = load_model(args.model, device)
    print(f"✅ Model loaded on {device}")
    
    # Load index
    print("Loading FAISS index...")
    index = faiss.read_index(args.index)
    with open(args.metadata, 'rb') as f:
        metadata = pickle.load(f)
    print(f"✅ Index loaded ({index.ntotal:,} vectors, nprobe={index.nprobe})")
    
    # Load reference
    print("Loading reference genome...")
    genome_dict = {}
    chr_lengths = {}
    for record in SeqIO.parse(args.reference, "fasta"):
        genome_dict[record.id] = str(record.seq)
        chr_lengths[record.id] = len(record.seq)
    print(f"✅ Loaded {len(genome_dict)} chromosomes")
    
    # Create pipeline
    print("\nInitializing pipeline...")
    seeder = AdaptiveSeeder(
        model=model,
        index=index,
        metadata=metadata,
        min_seeds=5,
        max_seeds=16,
        window_size=512,
        top_k=32
    )
    
    aligner = FastAligner(genome_dict, mode='semi-global')
    extend = ExtendPhase(aligner)
    print("✅ Pipeline ready")
    
    # Process reads
    print("\nProcessing reads...")
    reads = list(SeqIO.parse(args.reads, "fastq"))
    print(f"Found {len(reads)} reads")
    print()
    
    all_alignments = []
    start_time = time.time()
    
    for i, record in enumerate(reads):
        if (i + 1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            print(f"  Processed {i+1}/{len(reads)} reads ({rate:.2f} reads/sec)...")
        
        read_id = record.id
        read_seq = str(record.seq)
        
        try:
            # Seeding
            candidates = seeder.align_read(read_seq, read_id=read_id, return_top_k=5)
            
            if candidates is None:
                # Unmapped
                all_alignments.append({
                    'read_id': read_id,
                    'read_seq': read_seq,
                    'results': []
                })
                continue
            
            # EXTEND phase (with secondary alignments if requested)
            result = extend.extend_and_score(read_seq, candidates, return_all=args.secondary)
            
            if result is None:
                # Failed alignment
                all_alignments.append({
                    'read_id': read_id,
                    'read_seq': read_seq,
                    'results': []
                })
                continue
            
            # Prepare results for SAM output
            results_list = [result]
            
            # Add secondaries if present
            if args.secondary and 'secondary_alignments' in result:
                results_list.extend(result['secondary_alignments'])
            
            all_alignments.append({
                'read_id': read_id,
                'read_seq': read_seq,
                'results': results_list
            })
            
        except Exception as e:
            print(f"  ⚠️  Error processing {read_id}: {e}")
            all_alignments.append({
                'read_id': read_id,
                'read_seq': read_seq,
                'results': []
            })
    
    elapsed = time.time() - start_time
    
    # Write SAM output
    print(f"\nWriting SAM output to {args.output}...")
    write_sam_with_secondaries(args.output, all_alignments, chr_lengths)
    print("✅ SAM file written")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    mapped = sum(1 for a in all_alignments if len(a['results']) > 0)
    unmapped = len(all_alignments) - mapped
    secondaries = sum(len(a['results']) - 1 for a in all_alignments if len(a['results']) > 1)
    
    print(f"Total reads:       {len(reads)}")
    print(f"Mapped:            {mapped} ({mapped/len(reads)*100:.1f}%)")
    print(f"Unmapped:          {unmapped} ({unmapped/len(reads)*100:.1f}%)")
    if args.secondary:
        print(f"Secondary alns:    {secondaries}")
    print()
    print(f"Time:              {elapsed:.1f}s")
    print(f"Speed:             {len(reads)/elapsed:.2f} reads/sec")
    print()
    print(f"Output:            {args.output}")
    print("✅ Complete!")


if __name__ == '__main__':
    main()
