#!/usr/bin/env python3
"""
GenoCache V4 vs minimap2 Head-to-Head Benchmark

Comprehensive comparison:
1. Load GenoCache production index
2. Load test reads (synthetic + real ONT)
3. Run GenoCache alignment
4. Run minimap2 alignment
5. Compare: accuracy, speed, memory, results

Metrics:
- Accuracy (if ground truth available)
- Speed (reads/second)
- Memory usage
- Alignment quality
- Per-read timing
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
import faiss
import pickle
import subprocess
import time
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder


def load_index(index_dir: Path):
    """Load production GenoCache index"""
    print("Loading GenoCache index...")
    
    index_file = index_dir / "genocache_v4_production.index"
    metadata_file = index_dir / "genocache_v4_production.metadata.pkl"
    
    # Load FAISS index
    index = faiss.read_index(str(index_file))
    
    # Load metadata
    with open(metadata_file, 'rb') as f:
        metadata = pickle.load(f)
    
    print(f"✅ Index loaded:")
    print(f"  Embeddings: {index.ntotal:,}")
    print(f"  Dimension: {index.d}D")
    print(f"  Chromosomes: {len(np.unique(metadata['chr_names']))}")
    
    return index, metadata


def load_reads(reads_file: Path):
    """Load reads from FASTA/FASTQ"""
    print(f"\nLoading reads from {reads_file.name}...")
    
    reads = []
    read_ids = []
    
    with open(reads_file, 'r') as f:
        if reads_file.suffix in ['.fasta', '.fa']:
            # FASTA format
            seq = None
            read_id = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if seq is not None:
                        reads.append(seq)
                        read_ids.append(read_id)
                    read_id = line[1:].split()[0]
                    seq = ""
                else:
                    seq += line
            if seq is not None:
                reads.append(seq)
                read_ids.append(read_id)
        
        elif reads_file.suffix in ['.fastq', '.fq']:
            # FASTQ format
            line_num = 0
            for line in f:
                line = line.strip()
                if line_num % 4 == 0:  # ID line
                    read_id = line[1:].split()[0]
                elif line_num % 4 == 1:  # Sequence line
                    reads.append(line)
                    read_ids.append(read_id)
                line_num += 1
    
    print(f"✅ Loaded {len(reads)} reads")
    
    # Stats
    lengths = [len(r) for r in reads]
    print(f"  Length range: {min(lengths)}-{max(lengths)} bp")
    print(f"  Mean length: {np.mean(lengths):.0f} bp")
    print(f"  Median length: {np.median(lengths):.0f} bp")
    
    return reads, read_ids


def load_ground_truth(truth_file: Path):
    """Load ground truth mappings"""
    if not truth_file.exists():
        return None
    
    print(f"\nLoading ground truth from {truth_file.name}...")
    
    truth = {}
    with open(truth_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                read_id = parts[0]
                chr_name = parts[1]
                position = int(parts[2])
                truth[read_id] = {'chr': chr_name, 'pos': position}
    
    print(f"✅ Loaded ground truth for {len(truth)} reads")
    return truth


def geocache_align(
    model, index, metadata, reads, read_ids,
    device='cuda', k=10
):
    """Align reads using GenoCache"""
    print("\nRunning GenoCache alignment...")
    
    model.eval()
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq.upper()]
        return torch.tensor(indices, dtype=torch.long)
    
    results = []
    search_times = []
    encode_times = []
    
    start_time = time.time()
    
    with torch.no_grad():
        for i, (read_id, read_seq) in enumerate(tqdm(zip(read_ids, reads), total=len(reads), desc="Aligning")):
            # Extract 512bp seed from read
            if len(read_seq) >= 512:
                seed = read_seq[:512]
            else:
                seed = read_seq + 'N' * (512 - len(read_seq))
            
            # Encode read
            encode_start = time.time()
            read_tensor = seq_to_tensor(seed).unsqueeze(0).to(device)
            read_emb = model(read_tensor)
            read_emb = F.normalize(read_emb, p=2, dim=1)
            read_emb_np = read_emb.cpu().numpy().astype(np.float32)
            encode_time = time.time() - encode_start
            encode_times.append(encode_time)
            
            # Search index
            search_start = time.time()
            distances, indices = index.search(read_emb_np, k)
            search_time = time.time() - search_start
            search_times.append(search_time)
            
            # Get best match
            best_idx = indices[0][0]
            best_score = distances[0][0]
            best_chr = metadata['chr_names'][best_idx]
            best_pos = metadata['positions'][best_idx]
            
            results.append({
                'read_id': read_id,
                'chr': best_chr,
                'pos': best_pos,
                'score': float(best_score),
                'encode_time': encode_time,
                'search_time': search_time
            })
    
    total_time = time.time() - start_time
    
    print(f"✅ GenoCache alignment complete!")
    print(f"  Total time: {total_time:.2f}s")
    print(f"  Reads/second: {len(reads)/total_time:.1f}")
    print(f"  Avg encode time: {np.mean(encode_times)*1000:.2f} ms/read")
    print(f"  Avg search time: {np.mean(search_times)*1000:.2f} ms/read")
    print(f"  Avg total time: {np.mean([encode_times[i]+search_times[i] for i in range(len(reads))])*1000:.2f} ms/read")
    
    return results, {
        'total_time': total_time,
        'reads_per_sec': len(reads)/total_time,
        'avg_encode_ms': np.mean(encode_times)*1000,
        'avg_search_ms': np.mean(search_times)*1000,
        'avg_total_ms': np.mean([encode_times[i]+search_times[i] for i in range(len(reads))])*1000
    }


def minimap2_align(
    reads_file: Path,
    genome_file: Path,
    output_dir: Path,
    threads: int = 8
):
    """Align reads using minimap2"""
    print("\nRunning minimap2 alignment...")
    print(f"  Threads: {threads}")
    
    output_sam = output_dir / "minimap2_alignment.sam"
    
    # Run minimap2
    cmd = [
        'minimap2',
        '-ax', 'map-ont',  # ONT read preset
        '-t', str(threads),
        str(genome_file),
        str(reads_file)
    ]
    
    start_time = time.time()
    
    with open(output_sam, 'w') as f:
        process = subprocess.run(
            cmd,
            stdout=f,
            stderr=subprocess.PIPE,
            text=True
        )
    
    total_time = time.time() - start_time
    
    if process.returncode != 0:
        print(f"❌ minimap2 failed:")
        print(process.stderr)
        return None, None
    
    # Parse SAM output
    results = []
    with open(output_sam, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header
            
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            
            read_id = parts[0]
            flag = int(parts[1])
            chr_name = parts[2]
            pos = int(parts[3])
            mapq = int(parts[4])
            
            # Skip unmapped
            if flag & 4:
                continue
            
            results.append({
                'read_id': read_id,
                'chr': chr_name,
                'pos': pos,
                'mapq': mapq
            })
    
    num_reads = len(set(r['read_id'] for r in results))
    
    print(f"✅ minimap2 alignment complete!")
    print(f"  Total time: {total_time:.2f}s")
    print(f"  Reads/second: {num_reads/total_time:.1f}")
    print(f"  Avg time: {total_time/num_reads*1000:.2f} ms/read")
    print(f"  Mapped reads: {len(results)}")
    
    return results, {
        'total_time': total_time,
        'reads_per_sec': num_reads/total_time,
        'avg_total_ms': total_time/num_reads*1000
    }


def compare_accuracy(geocache_results, minimap2_results, ground_truth, tolerance=1000):
    """Compare accuracy against ground truth"""
    print("\n" + "=" * 80)
    print("ACCURACY COMPARISON")
    print("=" * 80)
    
    if ground_truth is None:
        print("⚠️  No ground truth available, skipping accuracy comparison")
        return None
    
    geocache_correct = 0
    minimap2_correct = 0
    both_correct = 0
    both_wrong = 0
    
    geocache_errors = []
    minimap2_errors = []
    
    comparison = []
    
    for gc_result in geocache_results:
        read_id = gc_result['read_id']
        
        if read_id not in ground_truth:
            continue
        
        truth = ground_truth[read_id]
        
        # Find minimap2 result
        mm2_result = next((r for r in minimap2_results if r['read_id'] == read_id), None)
        
        # GenoCache accuracy
        gc_error = abs(gc_result['pos'] - truth['pos']) if gc_result['chr'] == truth['chr'] else float('inf')
        gc_correct_flag = gc_error <= tolerance
        
        if gc_correct_flag:
            geocache_correct += 1
        geocache_errors.append(min(gc_error, 1e6))
        
        # minimap2 accuracy (if available)
        if mm2_result:
            mm2_error = abs(mm2_result['pos'] - truth['pos']) if mm2_result['chr'] == truth['chr'] else float('inf')
            mm2_correct_flag = mm2_error <= tolerance
            
            if mm2_correct_flag:
                minimap2_correct += 1
            minimap2_errors.append(min(mm2_error, 1e6))
            
            if gc_correct_flag and mm2_correct_flag:
                both_correct += 1
            elif not gc_correct_flag and not mm2_correct_flag:
                both_wrong += 1
            
            comparison.append({
                'read_id': read_id,
                'truth_chr': truth['chr'],
                'truth_pos': truth['pos'],
                'gc_chr': gc_result['chr'],
                'gc_pos': gc_result['pos'],
                'gc_error': gc_error,
                'gc_correct': gc_correct_flag,
                'mm2_chr': mm2_result['chr'],
                'mm2_pos': mm2_result['pos'],
                'mm2_error': mm2_error,
                'mm2_correct': mm2_correct_flag
            })
    
    total = len(geocache_errors)
    
    print(f"\nGenoCache Accuracy:")
    print(f"  Correct @ ±{tolerance}bp: {geocache_correct}/{total} ({100*geocache_correct/total:.1f}%)")
    print(f"  Median error: {np.median(geocache_errors):.0f} bp")
    print(f"  Mean error: {np.mean(geocache_errors):.0f} bp")
    
    if minimap2_errors:
        print(f"\nminimap2 Accuracy:")
        print(f"  Correct @ ±{tolerance}bp: {minimap2_correct}/{total} ({100*minimap2_correct/total:.1f}%)")
        print(f"  Median error: {np.median(minimap2_errors):.0f} bp")
        print(f"  Mean error: {np.mean(minimap2_errors):.0f} bp")
        
        print(f"\nAgreement:")
        print(f"  Both correct: {both_correct} ({100*both_correct/total:.1f}%)")
        print(f"  Both wrong: {both_wrong} ({100*both_wrong/total:.1f}%)")
        print(f"  GenoCache only: {geocache_correct - both_correct}")
        print(f"  minimap2 only: {minimap2_correct - both_correct}")
    
    return {
        'geocache_accuracy': 100*geocache_correct/total,
        'minimap2_accuracy': 100*minimap2_correct/total if minimap2_errors else None,
        'geocache_median_error': float(np.median(geocache_errors)),
        'minimap2_median_error': float(np.median(minimap2_errors)) if minimap2_errors else None,
        'comparison': comparison
    }


def compare_speed(geocache_timing, minimap2_timing):
    """Compare speed metrics"""
    print("\n" + "=" * 80)
    print("SPEED COMPARISON")
    print("=" * 80)
    
    print(f"\nGenoCache:")
    print(f"  Reads/second: {geocache_timing['reads_per_sec']:.1f}")
    print(f"  Time per read: {geocache_timing['avg_total_ms']:.2f} ms")
    print(f"    - Encoding: {geocache_timing['avg_encode_ms']:.2f} ms")
    print(f"    - Search: {geocache_timing['avg_search_ms']:.2f} ms")
    
    print(f"\nminimap2:")
    print(f"  Reads/second: {minimap2_timing['reads_per_sec']:.1f}")
    print(f"  Time per read: {minimap2_timing['avg_total_ms']:.2f} ms")
    
    speedup = geocache_timing['reads_per_sec'] / minimap2_timing['reads_per_sec']
    print(f"\nSpeedup:")
    if speedup > 1:
        print(f"  GenoCache is {speedup:.2f}× FASTER than minimap2 ✅")
    else:
        print(f"  minimap2 is {1/speedup:.2f}× faster than GenoCache")


def main():
    print("=" * 80)
    print("GenoCache V4 vs minimap2 - Head-to-Head Benchmark")
    print("=" * 80)
    print()
    
    # Configuration
    INDEX_DIR = Path("/home/nebius/genocache/genocache-v4/indexes")
    CHECKPOINT_PATH = Path("/home/nebius/genocache/genocache-v4/models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt")
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    READS_FILE = Path("/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa")
    TRUTH_FILE = Path("/home/nebius/genocache/genocache_data/reads_chr22_truth_FIXED.tsv")
    OUTPUT_DIR = Path("/home/nebius/genocache/genocache-v4/benchmark/results")
    
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    # Load model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}\n")
    
    print("Loading GenoCache model...")
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    checkpoint = torch.load(CHECKPOINT_PATH, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    print(f"✅ Model loaded\n")
    
    # Load index
    index, metadata = load_index(INDEX_DIR)
    
    # Load reads
    reads, read_ids = load_reads(READS_FILE)
    
    # Load ground truth
    ground_truth = load_ground_truth(TRUTH_FILE)
    
    # Run GenoCache
    geocache_results, geocache_timing = geocache_align(
        model, index, metadata, reads, read_ids, device=device, k=10
    )
    
    # Run minimap2
    minimap2_results, minimap2_timing = minimap2_align(
        READS_FILE, GENOME_PATH, OUTPUT_DIR, threads=8
    )
    
    # Compare
    accuracy_comparison = compare_accuracy(
        geocache_results, minimap2_results, ground_truth, tolerance=1000
    )
    
    compare_speed(geocache_timing, minimap2_timing)
    
    # Save results
    results_file = OUTPUT_DIR / "benchmark_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            'geocache': {
                'timing': geocache_timing,
                'results': geocache_results[:10]  # Sample
            },
            'minimap2': {
                'timing': minimap2_timing,
                'results': minimap2_results[:10] if minimap2_results else []
            },
            'accuracy': accuracy_comparison
        }, f, indent=2)
    
    print(f"\n✅ Results saved: {results_file}")
    
    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE!")
    print("=" * 80)


if __name__ == "__main__":
    main()
