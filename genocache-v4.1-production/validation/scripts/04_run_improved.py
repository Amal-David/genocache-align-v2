#!/usr/bin/env python3
"""
Run Improved Test (Current Code - With Improvements)

Tests current adaptive_seeding.py (with chain scoring + rescue improvements)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import torch
import faiss
import pickle
from Bio import SeqIO
import time

# Import modules
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "genocache_core"))
from encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner

def main():
    print("="*80)
    print("IMPROVED TEST (Current Code - With Improvements)")
    print("="*80)
    print()
    
    print("Using IMPROVED code:")
    print("  ✅ Chain scoring: Count anchors (fair)")
    print("  ✅ Rescue logic: 3 conditions (smart)")
    print("  ✅ nprobe: Preserved optimal value")
    print()
    
    # STEP 1: Load model and index
    print("Step 1: Loading model and index...")
    base_dir = Path(__file__).parent.parent.parent
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    checkpoint = torch.load(base_dir / "models" / "genocache_model.pt", map_location=device)
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model.load_state_dict(checkpoint)
    model.eval()
    print("  ✅ Model loaded")
    
    # Load index
    index = faiss.read_index(str(base_dir / "indexes" / "genocache_v4_production.index"))
    with open(base_dir / "indexes" / "genocache_v4_production.metadata.pkl", 'rb') as f:
        metadata = pickle.load(f)
    print(f"  ✅ Index loaded ({index.ntotal:,} vectors, nprobe={index.nprobe})")
    
    # Load genome
    print("  Loading reference...")
    genome_dict = {}
    ref_path = Path("/home/nebius/genocache/GRCh38.fa")
    for record in SeqIO.parse(ref_path, "fasta"):
        genome_dict[record.id] = str(record.seq)
    print(f"  ✅ Loaded {len(genome_dict)} chromosomes")
    
    # STEP 2: Create pipeline
    print("\nStep 2: Creating pipeline (IMPROVED VERSION)...")
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
    print("  ✅ Pipeline ready")
    
    # STEP 3: Load test data
    print("\nStep 3: Loading test data...")
    reads_file = Path(__file__).parent.parent / "data" / "test_reads_100.fastq"
    truth_file = Path(__file__).parent.parent / "data" / "test_reads_100_ground_truth.txt"
    
    reads = {}
    for record in SeqIO.parse(reads_file, "fastq"):
        reads[record.id] = str(record.seq)
    print(f"  ✅ Loaded {len(reads)} reads")
    
    # Load ground truth
    ground_truth = {}
    with open(truth_file) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            ground_truth[parts[0]] = {
                'chr': parts[1],
                'start': int(parts[2]),
                'end': int(parts[3])
            }
    print(f"  ✅ Loaded ground truth for {len(ground_truth)} reads")
    
    # STEP 4: Run test
    print("\nStep 4: Running IMPROVED test...")
    print("  (Using: count anchors + enhanced rescue)")
    print()
    
    results = []
    correct = 0
    wrong = 0
    unmapped = 0
    
    start_time = time.time()
    
    for i, (read_id, read_seq) in enumerate(reads.items()):
        if (i + 1) % 20 == 0:
            print(f"  Processed {i+1}/{len(reads)} reads...")
        
        try:
            # Run pipeline
            candidates = seeder.align_read(read_seq, read_id=read_id, return_top_k=5)
            
            if candidates is None:
                unmapped += 1
                results.append({
                    'read_id': read_id,
                    'predicted_chr': 'unmapped',
                    'true_chr': ground_truth[read_id]['chr'],
                    'correct': False,
                    'score': 0
                })
                continue
            
            # EXTEND phase
            result = extend.extend_and_score(read_seq, candidates)
            
            if result is None:
                unmapped += 1
                results.append({
                    'read_id': read_id,
                    'predicted_chr': 'unmapped',
                    'true_chr': ground_truth[read_id]['chr'],
                    'correct': False,
                    'score': 0
                })
                continue
            
            # Check correctness
            predicted_chr = result['chr']
            true_chr = ground_truth[read_id]['chr']
            is_correct = (predicted_chr == true_chr)
            
            if is_correct:
                correct += 1
            else:
                wrong += 1
            
            results.append({
                'read_id': read_id,
                'predicted_chr': predicted_chr,
                'true_chr': true_chr,
                'correct': is_correct,
                'score': result['alignment_score']
            })
            
        except Exception as e:
            print(f"  ⚠️  Error on {read_id}: {e}")
            unmapped += 1
            results.append({
                'read_id': read_id,
                'predicted_chr': 'error',
                'true_chr': ground_truth[read_id]['chr'],
                'correct': False,
                'score': 0
            })
    
    elapsed = time.time() - start_time
    
    # STEP 5: Save results
    print(f"\nStep 5: Saving results...")
    output_file = Path(__file__).parent.parent / "results" / "improved_results.txt"
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("read_id\tpredicted_chr\ttrue_chr\tcorrect\tscore\n")
        for r in results:
            f.write(f"{r['read_id']}\t{r['predicted_chr']}\t{r['true_chr']}\t{r['correct']}\t{r['score']}\n")
    
    print(f"  ✅ Results saved to: {output_file}")
    
    # STEP 6: Summary
    print("\n" + "="*80)
    print("IMPROVED RESULTS")
    print("="*80)
    
    total = len(reads)
    accuracy = (correct / total) * 100 if total > 0 else 0
    
    print(f"\nTotal reads:   {total}")
    print(f"Correct:       {correct} ({accuracy:.1f}%)")
    print(f"Wrong:         {wrong} ({wrong/total*100 if total > 0 else 0:.1f}%)")
    print(f"Unmapped:      {unmapped} ({unmapped/total*100 if total > 0 else 0:.1f}%)")
    print()
    print(f"Speed:         {elapsed:.1f}s total ({total/elapsed:.2f} reads/sec)")
    print()
    print(f"Accuracy:      {accuracy:.1f}%")


if __name__ == '__main__':
    main()
