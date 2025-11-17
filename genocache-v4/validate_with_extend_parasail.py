#!/usr/bin/env python3
"""
Validate EXTEND Phase Fix with Parasail

This script validates that the EXTEND phase fixes our 37% chromosome accuracy bug.

OLD Method (WRONG):
- Pick candidate by seed count
- Result: 37.5% chromosome accuracy (3/8 correct)

NEW Method (CORRECT - with EXTEND):
- Get top-k candidates from seeding
- Align to EACH candidate using parasail
- Pick best by ALIGNMENT SCORE
- Result: Expected 95%+ chromosome accuracy

This is the critical missing piece from NeuralAligner!
"""

import sys
import torch
import faiss
import parasail
import numpy as np
from pathlib import Path
from Bio import SeqIO

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from adaptive_seeding import AdaptiveSeeder
from fast_alignment import FastAligner
from extend_phase import ExtendPhase
from models.encoder import GenoCacheEncoder

def load_model(checkpoint_path):
    """Load GenoCache model"""
    print("Loading GenoCache model...")
    model = GenoCacheEncoder(
        emb_dim=128,  # Actual trained model uses 128D
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    )
    
    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    elif 'state_dict' in checkpoint:
        model.load_state_dict(checkpoint['state_dict'])
    else:
        model.load_state_dict(checkpoint)
    
    model.eval()
    print(f"✅ Model loaded from {checkpoint_path}")
    return model

def load_index_and_metadata(index_path, metadata_path):
    """Load FAISS index and metadata"""
    print("Loading FAISS index...")
    index = faiss.read_index(str(index_path))
    print(f"✅ Index loaded: {index.ntotal:,} vectors")
    
    print("Loading metadata...")
    import pickle
    with open(metadata_path, 'rb') as f:
        metadata = pickle.load(f)
    print(f"✅ Metadata loaded: {len(metadata['chr_names']):,} entries")
    
    return index, metadata

def load_reference_genome(fasta_path, chromosomes=None):
    """Load reference genome"""
    print("Loading reference genome...")
    genome_dict = {}
    
    # Load all main chromosomes (chr1-22, X, Y)
    # Don't restrict - load everything for validation
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Load all chromosomes (main + alternate)
        genome_dict[record.id] = str(record.seq)
        if len(genome_dict) <= 30:  # Print first 30
            print(f"  Loaded {record.id}: {len(record.seq):,} bp")
    
    if len(genome_dict) > 30:
        print(f"  ... and {len(genome_dict) - 30} more")
    
    print(f"✅ Loaded {len(genome_dict)} chromosomes/contigs")
    return genome_dict

def parse_sam_alignments(sam_file):
    """Parse SAM file to get read sequences and expected positions"""
    reads = {}
    
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            read_id = fields[0]
            flag = int(fields[1])
            ref_chr = fields[2]
            pos = int(fields[3])
            seq = fields[9]
            
            # Skip unmapped
            if flag & 4:
                continue
            
            reads[read_id] = {
                'seq': seq,
                'expected_chr': ref_chr,
                'expected_pos': pos
            }
    
    return reads

def align_with_extend_phase(read_id, read_seq, seeder, aligner, extend_phase, 
                           top_k=5):
    """
    Align read using EXTEND phase
    
    Args:
        read_id: Read identifier
        read_seq: Read sequence
        seeder: AdaptiveSeeder instance
        aligner: FastAligner instance
        extend_phase: ExtendPhase instance
        top_k: Number of candidates to test
    
    Returns:
        Best alignment result
    """
    print(f"\n{'='*80}")
    print(f"Processing {read_id}")
    print(f"{'='*80}")
    
    # Step 1: Adaptive seeding (get top-k candidates)
    print(f"Step 1: Adaptive seeding...")
    candidates = seeder.align_read(read_seq, read_id, return_top_k=top_k)
    
    if not candidates:
        print(f"❌ Seeding failed - no candidates found")
        return None
    
    print(f"  → Found {len(candidates)} candidate chains")
    
    # Show candidates
    for i, candidate in enumerate(candidates):
        print(f"    Candidate {i+1}: {candidate['chr']} "
              f"({candidate['num_seeds']} seeds, score={candidate['score']:.2f})")
    
    # Step 2: EXTEND phase - align to each candidate, pick best by score
    print(f"\nStep 2: EXTEND phase (align to each candidate)...")
    result = extend_phase.extend_and_score(read_seq, candidates)
    
    if not result:
        print(f"❌ EXTEND phase failed - all alignments below threshold")
        return None
    
    # Show results
    print(f"\n✅ EXTEND phase complete:")
    print(f"  Best chromosome: {result['chr']}")
    print(f"  Alignment score: {result['alignment_score']}")
    print(f"  Position: {result['start']:,} - {result['end']:,}")
    print(f"  CIGAR: {result['cigar']}")
    print(f"  Status: {result['status']}")
    print(f"  Candidates tested: {result['candidates_tested']}")
    print(f"  Candidates aligned: {result['candidates_aligned']}")
    
    return result

def main():
    print("="*80)
    print("VALIDATING EXTEND PHASE FIX WITH PARASAIL")
    print("="*80)
    print()
    
    # Paths
    base_dir = Path("/home/nebius/genocache/genocache-v4")
    model_path = base_dir / "models/checkpoints/fullgenome_best_sep12.3547_epoch8.pt"
    index_path = base_dir / "indexes/genocache_v4_production.index"
    metadata_path = base_dir / "indexes/genocache_v4_production.metadata.pkl"
    genome_path = Path("/home/nebius/genocache/GRCh38.fa")
    old_sam = base_dir / "test_complete_10reads.sam"
    minimap2_sam = base_dir / "minimap2_same_10reads.sam"
    
    # Check files
    if not model_path.exists():
        print(f"❌ Model not found: {model_path}")
        return
    if not index_path.exists():
        print(f"❌ Index not found: {index_path}")
        return
    if not genome_path.exists():
        print(f"❌ Reference not found: {genome_path}")
        return
    if not old_sam.exists():
        print(f"❌ OLD SAM not found: {old_sam}")
        return
    if not minimap2_sam.exists():
        print(f"❌ minimap2 SAM not found: {minimap2_sam}")
        return
    
    # Load components
    model = load_model(model_path)
    index, metadata = load_index_and_metadata(index_path, metadata_path)
    
    # Load reference (only needed chromosomes)
    genome_dict = load_reference_genome(genome_path)
    
    # Initialize components
    print("\nInitializing components...")
    seeder = AdaptiveSeeder(
        model=model,
        index=index,
        metadata=metadata,
        min_seeds=5,
        max_seeds=16,
        window_size=512,
        top_k=32
    )
    print("✅ Adaptive seeder initialized")
    
    aligner = FastAligner(genome_dict=genome_dict, mode='semi-global')
    print("✅ Fast aligner initialized")
    
    extend_phase = ExtendPhase(
        aligner=aligner,
        min_score_threshold=100,
        score_ratio_threshold=1.5
    )
    print("✅ EXTEND phase initialized")
    
    # Parse OLD and minimap2 results
    print("\nParsing OLD results and ground truth...")
    old_results = parse_sam_alignments(old_sam)
    minimap2_results = parse_sam_alignments(minimap2_sam)
    
    print(f"OLD results: {len(old_results)} reads")
    print(f"minimap2 (ground truth): {len(minimap2_results)} reads")
    
    # Get common reads
    common_reads = set(old_results.keys()) & set(minimap2_results.keys())
    print(f"Common reads: {len(common_reads)}")
    
    # Test EXTEND phase on each read
    print("\n" + "="*80)
    print("TESTING EXTEND PHASE")
    print("="*80)
    
    new_results = {}
    
    for read_id in sorted(common_reads):
        read_info = old_results[read_id]
        read_seq = read_info['seq']
        
        # Align with EXTEND phase
        result = align_with_extend_phase(
            read_id, read_seq, seeder, aligner, extend_phase, top_k=5
        )
        
        if result:
            new_results[read_id] = result
    
    # Compare results
    print("\n" + "="*80)
    print("COMPARISON: OLD vs NEW vs GROUND TRUTH")
    print("="*80)
    print()
    
    old_chr_correct = 0
    new_chr_correct = 0
    total = len(common_reads)
    
    for read_id in sorted(common_reads):
        ground_truth_chr = minimap2_results[read_id]['expected_chr']
        old_chr = old_results[read_id]['expected_chr']
        
        if read_id in new_results:
            new_chr = new_results[read_id]['chr']
        else:
            new_chr = "UNMAPPED"
        
        old_correct = (old_chr == ground_truth_chr)
        new_correct = (new_chr == ground_truth_chr)
        
        if old_correct:
            old_chr_correct += 1
        if new_correct:
            new_chr_correct += 1
        
        print(f"{read_id}:")
        print(f"  Ground truth: {ground_truth_chr}")
        print(f"  OLD method:   {old_chr} {'✅' if old_correct else '❌'}")
        print(f"  NEW (EXTEND): {new_chr} {'✅' if new_correct else '❌'}")
        
        if not old_correct and new_correct:
            print(f"  → ⭐ FIXED by EXTEND phase!")
        elif old_correct and not new_correct:
            print(f"  → ⚠️  REGRESSED!")
        print()
    
    # Final summary
    print("="*80)
    print("FINAL RESULTS")
    print("="*80)
    
    old_accuracy = old_chr_correct / total * 100
    new_accuracy = new_chr_correct / total * 100
    improvement = new_accuracy - old_accuracy
    
    print(f"Total reads tested: {total}")
    print()
    print(f"OLD Method (seed count):")
    print(f"  Correct: {old_chr_correct}/{total} ({old_accuracy:.1f}%)")
    print()
    print(f"NEW Method (EXTEND + alignment scores):")
    print(f"  Correct: {new_chr_correct}/{total} ({new_accuracy:.1f}%)")
    print()
    print(f"Improvement: {improvement:+.1f}% ({'+' if improvement > 0 else ''}{new_chr_correct - old_chr_correct} reads)")
    print()
    
    if improvement > 30:
        print("🎉 MAJOR IMPROVEMENT! EXTEND phase works!")
    elif improvement > 10:
        print("✅ Good improvement! EXTEND phase helps!")
    elif improvement > 0:
        print("✓ Small improvement")
    else:
        print("⚠️  No improvement - need to debug")
    
    print()
    print("="*80)

if __name__ == '__main__':
    main()
