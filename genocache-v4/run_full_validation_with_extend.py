#!/usr/bin/env python3
"""
Full Pipeline Validation with EXTEND Phase

This script:
1. Loads the trained model and FAISS index
2. Runs complete pipeline with EXTEND phase on test reads
3. Generates proper SAM output with CIGAR, MAPQ, etc.
4. Compares with OLD method and minimap2

Expected: 37% → 95%+ chromosome accuracy improvement
"""

import sys
import time
from pathlib import Path
import torch
import faiss
import pickle
import numpy as np

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from model import GenoCache
from adaptive_seeding import AdaptiveSeeder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner

def load_fasta(fasta_file):
    """Load reads from FASTA file"""
    reads = []
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    reads.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            reads.append((current_id, ''.join(current_seq)))
    
    return reads

def calculate_mapq(alignment_score, second_best_score=0):
    """Calculate MAPQ score from alignment scores"""
    if second_best_score == 0 or alignment_score <= second_best_score:
        # Primary alignment, high confidence
        return 60
    
    # MAPQ based on score difference
    score_diff = alignment_score - second_best_score
    if score_diff > 100:
        return 60
    elif score_diff > 50:
        return 40
    elif score_diff > 20:
        return 20
    else:
        return 10

def write_sam_header(sam_file, ref_file):
    """Write SAM header"""
    with open(sam_file, 'w') as f:
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        
        # Read reference to get chromosome lengths
        chr_lengths = {}
        current_chr = None
        current_len = 0
        
        with open(ref_file) as ref:
            for line in ref:
                if line.startswith('>'):
                    if current_chr:
                        chr_lengths[current_chr] = current_len
                    current_chr = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line.strip())
            
            if current_chr:
                chr_lengths[current_chr] = current_len
        
        # Write chromosome headers
        for chr_name, length in sorted(chr_lengths.items())[:25]:  # Main chromosomes
            f.write(f"@SQ\tSN:{chr_name}\tLN:{length}\n")
        
        f.write("@PG\tID:genocache\tPN:genocache\tVN:4.0\tCL:genocache_with_extend\n")

def write_sam_alignment(sam_file, read_id, read_seq, alignment_result):
    """Write alignment to SAM file"""
    with open(sam_file, 'a') as f:
        if alignment_result is None:
            # Unmapped
            f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
            return
        
        # Mapped alignment
        flag = 0  # Primary alignment
        chr_name = alignment_result['chr']
        pos = alignment_result['start'] + 1  # SAM is 1-based
        mapq = alignment_result.get('mapq', 60)
        cigar = alignment_result.get('cigar', f"{len(read_seq)}M")  # Default to match
        
        # Optional fields
        tags = []
        tags.append(f"AS:i:{alignment_result.get('alignment_score', 0)}")
        tags.append(f"XS:i:{alignment_result.get('num_seeds', 0)}")
        tags.append(f"XC:i:{alignment_result.get('candidates_tested', 0)}")
        
        tags_str = '\t'.join(tags)
        
        f.write(f"{read_id}\t{flag}\t{chr_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read_seq}\t*\t{tags_str}\n")

def main():
    print("=" * 80)
    print("FULL VALIDATION WITH EXTEND PHASE")
    print("=" * 80)
    print()
    
    # Configuration
    model_path = 'models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt'
    index_path = 'indexes/genocache_v4_production.index'
    metadata_path = 'indexes/genocache_v4_production.metadata.pkl'
    ref_path = '../GRCh38.fa'
    test_reads = 'test_10_reads_exact.fa'
    output_sam = 'test_with_extend_phase.sam'
    
    # Check files exist
    for path, name in [(model_path, 'Model'), (index_path, 'Index'), 
                       (metadata_path, 'Metadata'), (ref_path, 'Reference'),
                       (test_reads, 'Test reads')]:
        if not Path(path).exists():
            print(f"❌ Error: {name} not found: {path}")
            return 1
    
    print("Loading components...")
    print()
    
    # Load model
    print(f"1. Loading model from {model_path}...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"   Using device: {device}")
    
    checkpoint = torch.load(model_path, map_location=device, weights_only=True)
    
    # Get model config from checkpoint
    model_config = checkpoint.get('model_config', {
        'vocab_size': 5,
        'embed_dim': 128,
        'num_layers': 6,
        'max_len': 512
    })
    
    model = GenoCache(**model_config).to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    print(f"   ✅ Model loaded ({model_config.get('d_model', 128)}D embeddings)")
    print()
    
    # Load FAISS index
    print(f"2. Loading FAISS index from {index_path}...")
    index = faiss.read_index(index_path)
    print(f"   ✅ Index loaded ({index.ntotal:,} vectors)")
    print()
    
    # Load metadata
    print(f"3. Loading metadata from {metadata_path}...")
    with open(metadata_path, 'rb') as f:
        metadata = pickle.load(f)
    print(f"   ✅ Metadata loaded ({len(metadata):,} entries)")
    print()
    
    # Initialize components
    print("4. Initializing pipeline components...")
    
    # Adaptive seeder (returns top-k candidates)
    seeder = AdaptiveSeeder(
        model=model,
        index=index,
        metadata=metadata,
        min_seeds=5,
        max_seeds=16,
        window_size=512,
        top_k=32
    )
    print("   ✅ Adaptive seeder initialized (returns top-5 candidates)")
    
    # Fast aligner (for EXTEND phase)
    aligner = FastAligner(ref_path)
    print("   ✅ Fast aligner initialized (parasail)")
    
    # EXTEND phase (THE FIX!)
    extend = ExtendPhase(
        aligner=aligner,
        min_score_threshold=100,
        score_ratio_threshold=1.5
    )
    print("   ✅ EXTEND phase initialized (align to each candidate)")
    print()
    
    # Load test reads
    print(f"5. Loading test reads from {test_reads}...")
    reads = load_fasta(test_reads)
    print(f"   ✅ Loaded {len(reads)} reads")
    print()
    
    # Initialize SAM output
    print(f"6. Initializing SAM output: {output_sam}...")
    write_sam_header(output_sam, ref_path)
    print(f"   ✅ SAM header written")
    print()
    
    print("=" * 80)
    print("RUNNING PIPELINE WITH EXTEND PHASE")
    print("=" * 80)
    print()
    
    total_time = 0
    results = []
    
    for i, (read_id, read_seq) in enumerate(reads):
        print(f"Read {i+1}/{len(reads)}: {read_id}")
        print(f"  Length: {len(read_seq)} bp")
        
        start_time = time.time()
        
        try:
            # Step 1: Adaptive seeding (get top-k candidates)
            print(f"  [1/3] Adaptive seeding...")
            seed_start = time.time()
            
            candidates = seeder.get_top_k_candidates(read_seq, read_id, top_k=5)
            
            seed_time = time.time() - seed_start
            print(f"        → Found {len(candidates)} candidates ({seed_time:.3f}s)")
            
            if candidates:
                for j, cand in enumerate(candidates[:3]):
                    print(f"           {j+1}. {cand['chr']}:{cand['start']}-{cand['end']} "
                          f"({cand.get('num_seeds', 0)} seeds, score={cand.get('score', 0):.1f})")
            
            # Step 2: EXTEND phase (THE FIX!)
            print(f"  [2/3] EXTEND phase (align to each candidate)...")
            extend_start = time.time()
            
            if candidates:
                alignment_result = extend.extend_and_score(read_seq, candidates)
            else:
                alignment_result = None
            
            extend_time = time.time() - extend_start
            
            if alignment_result:
                print(f"        → Best: {alignment_result['chr']}:{alignment_result['start']}")
                print(f"           Alignment score: {alignment_result['alignment_score']}")
                print(f"           Status: {alignment_result['status']}")
                print(f"           Candidates tested: {alignment_result['candidates_tested']}")
                print(f"           Time: {extend_time:.3f}s")
                
                # Calculate MAPQ
                alignment_result['mapq'] = calculate_mapq(
                    alignment_result['alignment_score']
                )
            else:
                print(f"        → No alignment found")
            
            # Step 3: Write to SAM
            print(f"  [3/3] Writing to SAM...")
            write_sam_alignment(output_sam, read_id, read_seq, alignment_result)
            
            total_read_time = time.time() - start_time
            total_time += total_read_time
            
            print(f"  Total time: {total_read_time:.3f}s")
            print()
            
            results.append({
                'read_id': read_id,
                'mapped': alignment_result is not None,
                'chr': alignment_result['chr'] if alignment_result else None,
                'pos': alignment_result['start'] if alignment_result else None,
                'time': total_read_time
            })
            
        except Exception as e:
            print(f"  ❌ Error processing read: {e}")
            import traceback
            traceback.print_exc()
            write_sam_alignment(output_sam, read_id, read_seq, None)
            print()
    
    print("=" * 80)
    print("VALIDATION COMPLETE")
    print("=" * 80)
    print()
    
    # Summary statistics
    mapped_count = sum(1 for r in results if r['mapped'])
    avg_time = total_time / len(reads) if reads else 0
    
    print(f"Total reads: {len(reads)}")
    print(f"Mapped: {mapped_count}/{len(reads)} ({100*mapped_count/len(reads):.1f}%)")
    print(f"Unmapped: {len(reads)-mapped_count}")
    print(f"Average time per read: {avg_time:.3f}s")
    print(f"Total time: {total_time:.3f}s")
    print()
    
    print(f"✅ SAM output written to: {output_sam}")
    print()
    
    print("Next steps:")
    print(f"1. Compare with OLD method:")
    print(f"   python3 compare_chromosome_accuracy.py {output_sam} minimap2_same_10reads.sam")
    print()
    print(f"2. View SAM file:")
    print(f"   head -50 {output_sam}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
