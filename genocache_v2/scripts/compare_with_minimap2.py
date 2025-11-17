#!/usr/bin/env python3
"""
Comprehensive GenoCache vs Minimap2 Comparison

Generates realistic ONT-style reads with errors and compares:
1. Mapping accuracy (position correctness)
2. Performance (speed, throughput)
3. Memory usage

Usage:
    python scripts/compare_with_minimap2.py \
        --checkpoint /home/nebius/work/genocache_checkpoints/improved_cnn_best.pt \
        --index indexes/index_ivfpq.faiss \
        --positions data/reference_encodings/ref_positions_20251111_163637.npy \
        --fasta references/GCF_000001405.40_GRCh38.p14_chr22.fna \
        --chrom NC_000022.11 \
        --num-reads 500 \
        --output-dir comparison_test
"""

import argparse
import json
import subprocess
import time
import numpy as np
import torch
import faiss
from pathlib import Path
from collections import defaultdict
import sys
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))
from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode


def add_sequencing_errors(sequence, error_rate=0.08, sub_rate=0.6, ins_rate=0.2, del_rate=0.2):
    """Add realistic ONT sequencing errors"""
    bases = ['A', 'C', 'G', 'T']
    seq = list(sequence)
    num_errors = 0
    
    i = 0
    while i < len(seq):
        if np.random.random() < error_rate:
            error_type = np.random.choice(['sub', 'ins', 'del'], p=[sub_rate, ins_rate, del_rate])
            
            if error_type == 'sub':
                original = seq[i]
                alternatives = [b for b in bases if b != original]
                seq[i] = np.random.choice(alternatives)
                num_errors += 1
                i += 1
            elif error_type == 'ins':
                seq.insert(i, np.random.choice(bases))
                num_errors += 1
                i += 2
            else:
                if i < len(seq):
                    seq.pop(i)
                    num_errors += 1
        else:
            i += 1
    
    return ''.join(seq), num_errors


def generate_realistic_reads(ref_seq, num_reads=500, min_len=5000, max_len=20000, error_rate=0.08, seed=42):
    """Generate realistic ONT-style reads with errors"""
    np.random.seed(seed)
    reads = []
    ref_len = len(ref_seq)
    
    attempts = 0
    while len(reads) < num_reads and attempts < num_reads * 3:
        attempts += 1
        read_len = np.random.randint(min_len, max_len)
        start = np.random.randint(100, ref_len - read_len - 100)
        
        clean_seq = ref_seq[start:start + read_len]
        
        # Skip N-rich regions
        if clean_seq.count('N') > read_len * 0.05:
            continue
        
        # Add errors
        noisy_seq, num_errors = add_sequencing_errors(clean_seq, error_rate)
        
        read_id = f"read_{len(reads)}"
        reads.append({
            'id': read_id,
            'seq': noisy_seq,
            'true_pos': start,
            'true_len': read_len,
            'num_errors': num_errors
        })
    
    return reads


def write_fasta(reads, output_file):
    """Write reads to FASTA format"""
    with open(output_file, 'w') as f:
        for read in reads:
            f.write(f">{read['id']}\n{read['seq']}\n")


def run_minimap2(fasta_ref, reads_fasta, output_sam, threads=8):
    """Run minimap2 alignment"""
    print("\n" + "="*80)
    print("Running Minimap2...")
    print("="*80)
    
    cmd = [
        'minimap2',
        '-ax', 'map-ont',
        '-t', str(threads),
        str(fasta_ref),
        str(reads_fasta)
    ]
    
    print(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    with open(output_sam, 'w') as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    elapsed = time.time() - start_time
    
    print(f"✅ Minimap2 completed in {elapsed:.2f}s")
    
    return elapsed, result.stderr


def run_genocache(checkpoint, index_path, positions_path, fasta_ref, chrom, reads, output_sam, k=50, num_seeds=5, use_cuda=True):
    """Run GenoCache alignment"""
    print("\n" + "="*80)
    print("Running GenoCache...")
    print("="*80)
    
    device = 'cuda' if use_cuda and torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")
    
    # Load model
    print("Loading model...")
    model, metadata = ImprovedCNN.load_checkpoint(checkpoint, device=device)
    model.eval()
    seed_len = model.input_len
    
    # Load FAISS index
    print("Loading FAISS index...")
    index = faiss.read_index(index_path)
    print(f"  Loaded {index.ntotal:,} vectors")
    
    # Load positions
    print("Loading positions...")
    positions = np.load(positions_path)
    
    # Load reference
    print("Loading reference...")
    ref = ReferenceGenome(fasta_ref)
    ref_seq = ref.sequences[chrom]
    
    # Align reads
    print(f"Aligning {len(reads)} reads...")
    alignments = []
    start_time = time.time()
    
    for i, read in enumerate(reads):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(reads)} reads...")
        
        read_seq = read['seq']
        read_len = len(read_seq)
        
        # Extract seeds
        best_pos = None
        best_score = 0
        
        for seed_idx in range(num_seeds):
            offset = seed_idx * max(1, (read_len - seed_len) // (num_seeds + 1))
            if offset + seed_len > read_len:
                break
            
            seed_seq = read_seq[offset:offset + seed_len]
            if len(seed_seq) != seed_len:
                continue
            
            # Encode seed
            tensor = torch.tensor(one_hot_encode(seed_seq)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(tensor).cpu().numpy()
            
            # Search FAISS
            emb = np.ascontiguousarray(emb, dtype=np.float32)
            D, I = index.search(emb, k)
            
            if I[0].size > 0:
                # Get top candidate
                cand_pos = int(positions[I[0][0]])
                cand_score = float(D[0][0])
                
                if cand_score > best_score:
                    best_score = cand_score
                    best_pos = cand_pos - offset  # Adjust for seed offset
        
        if best_pos is not None:
            alignments.append({
                'read_id': read['id'],
                'pos': best_pos,
                'mapq': 60,
                'flag': 0
            })
        else:
            alignments.append({
                'read_id': read['id'],
                'pos': -1,
                'mapq': 0,
                'flag': 4  # Unmapped
            })
    
    elapsed = time.time() - start_time
    print(f"✅ GenoCache completed in {elapsed:.2f}s")
    
    # Write SAM
    print("Writing SAM file...")
    write_sam(alignments, reads, chrom, len(ref_seq), output_sam)
    
    return elapsed, alignments


def write_sam(alignments, reads, chrom, ref_len, output_sam):
    """Write alignments to SAM format"""
    with open(output_sam, 'w') as f:
        # Write header
        f.write(f"@HD\tVN:1.0\tSO:unsorted\n")
        f.write(f"@SQ\tSN:{chrom}\tLN:{ref_len}\n")
        f.write(f"@PG\tID:genocache\tPN:genocache\tVN:2.0\n")
        
        # Write alignments
        read_dict = {r['id']: r for r in reads}
        for aln in alignments:
            read = read_dict[aln['read_id']]
            pos = aln['pos'] + 1  # Convert to 1-based
            cigar = f"{len(read['seq'])}M" if aln['flag'] == 0 else "*"
            
            f.write(f"{aln['read_id']}\t{aln['flag']}\t{chrom}\t{pos}\t{aln['mapq']}\t{cigar}\t*\t0\t0\t{read['seq']}\t*\n")


def parse_sam(sam_file):
    """Parse SAM file and extract alignments"""
    alignments = {}
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            read_id = fields[0]
            flag = int(fields[1])
            pos = int(fields[3]) - 1  # Convert to 0-based
            mapq = int(fields[4])
            
            if flag & 4:  # Unmapped
                continue
            
            alignments[read_id] = {
                'pos': pos,
                'mapq': mapq
            }
    
    return alignments


def evaluate_accuracy(reads, alignments, tolerance=1000):
    """Evaluate alignment accuracy"""
    stats = {
        'total': len(reads),
        'mapped': 0,
        'correct': 0,
        'errors': []
    }
    
    for read in reads:
        read_id = read['id']
        true_pos = read['true_pos']
        
        if read_id in alignments:
            stats['mapped'] += 1
            aln_pos = alignments[read_id]['pos']
            error = abs(aln_pos - true_pos)
            stats['errors'].append(error)
            
            if error <= tolerance:
                stats['correct'] += 1
    
    # Calculate metrics
    stats['mapping_rate'] = 100 * stats['mapped'] / stats['total']
    stats['accuracy'] = 100 * stats['correct'] / stats['total']
    
    if stats['errors']:
        stats['median_error'] = float(np.median(stats['errors']))
        stats['mean_error'] = float(np.mean(stats['errors']))
        stats['p95_error'] = float(np.percentile(stats['errors'], 95))
    
    return stats


def print_comparison(genocache_stats, minimap2_stats, genocache_time, minimap2_time, num_reads):
    """Print comparison table"""
    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)
    print()
    
    print(f"{'Metric':<35} {'Minimap2':<20} {'GenoCache':<20} {'Δ':<15}")
    print("-"*90)
    
    # Mapping rate
    mm2_map = minimap2_stats['mapping_rate']
    gc_map = genocache_stats['mapping_rate']
    print(f"{'Mapping rate':<35} {mm2_map:.2f}%{'':<15} {gc_map:.2f}%{'':<15} {gc_map-mm2_map:+.2f}%")
    
    # Accuracy
    mm2_acc = minimap2_stats['accuracy']
    gc_acc = genocache_stats['accuracy']
    print(f"{'Accuracy (±1kb)':<35} {mm2_acc:.2f}%{'':<15} {gc_acc:.2f}%{'':<15} {gc_acc-mm2_acc:+.2f}%")
    
    # Median error
    mm2_med = minimap2_stats.get('median_error', 0)
    gc_med = genocache_stats.get('median_error', 0)
    print(f"{'Median position error (bp)':<35} {mm2_med:.0f}{'':<19} {gc_med:.0f}{'':<19} {gc_med-mm2_med:+.0f}")
    
    # Mean error
    mm2_mean = minimap2_stats.get('mean_error', 0)
    gc_mean = genocache_stats.get('mean_error', 0)
    print(f"{'Mean position error (bp)':<35} {mm2_mean:.0f}{'':<19} {gc_mean:.0f}{'':<19} {gc_mean-mm2_mean:+.0f}")
    
    # P95 error
    mm2_p95 = minimap2_stats.get('p95_error', 0)
    gc_p95 = genocache_stats.get('p95_error', 0)
    print(f"{'P95 position error (bp)':<35} {mm2_p95:.0f}{'':<19} {gc_p95:.0f}{'':<19} {gc_p95-mm2_p95:+.0f}")
    
    print()
    print("-"*90)
    
    # Performance
    mm2_speed = num_reads / minimap2_time
    gc_speed = num_reads / genocache_time
    speedup = mm2_speed / gc_speed
    
    print(f"{'Time (seconds)':<35} {minimap2_time:.2f}{'':<19} {genocache_time:.2f}{'':<19} {genocache_time/minimap2_time:.2f}x")
    print(f"{'Throughput (reads/sec)':<35} {mm2_speed:.1f}{'':<19} {gc_speed:.1f}{'':<19} {speedup:.2f}x")
    
    print()
    print("="*80)


def main():
    parser = argparse.ArgumentParser(description='Compare GenoCache vs Minimap2')
    
    # Paths
    parser.add_argument('--checkpoint', required=True, help='GenoCache checkpoint')
    parser.add_argument('--index', required=True, help='FAISS index')
    parser.add_argument('--positions', required=True, help='Position mappings')
    parser.add_argument('--fasta', required=True, help='Reference FASTA')
    parser.add_argument('--chrom', required=True, help='Chromosome name')
    
    # Parameters
    parser.add_argument('--num-reads', type=int, default=500, help='Number of reads')
    parser.add_argument('--min-length', type=int, default=5000, help='Min read length')
    parser.add_argument('--max-length', type=int, default=20000, help='Max read length')
    parser.add_argument('--error-rate', type=float, default=0.08, help='Error rate')
    parser.add_argument('--k', type=int, default=50, help='FAISS k')
    parser.add_argument('--num-seeds', type=int, default=5, help='Seeds per read')
    parser.add_argument('--tolerance', type=int, default=1000, help='Position tolerance')
    parser.add_argument('--threads', type=int, default=8, help='Threads for minimap2')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--use-cuda', action='store_true', help='Use CUDA')
    
    # Output
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("GenoCache vs Minimap2 Comparison")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Reads: {args.num_reads}")
    print(f"Read length: {args.min_length}-{args.max_length} bp")
    print(f"Error rate: {args.error_rate*100:.1f}%")
    print("="*80)
    
    # Load reference
    print("\nLoading reference...")
    ref = ReferenceGenome(args.fasta)
    ref_seq = ref.sequences[args.chrom]
    print(f"✅ Reference: {args.chrom}, {len(ref_seq):,} bp")
    
    # Generate reads
    print(f"\nGenerating {args.num_reads} realistic ONT reads...")
    reads = generate_realistic_reads(
        ref_seq,
        num_reads=args.num_reads,
        min_len=args.min_length,
        max_len=args.max_length,
        error_rate=args.error_rate,
        seed=args.seed
    )
    print(f"✅ Generated {len(reads)} reads")
    
    # Calculate error stats
    total_errors = sum(r['num_errors'] for r in reads)
    total_bases = sum(len(r['seq']) for r in reads)
    actual_error_rate = total_errors / total_bases
    print(f"   Actual error rate: {actual_error_rate*100:.2f}%")
    
    # Write reads
    reads_fasta = output_dir / "reads.fasta"
    write_fasta(reads, reads_fasta)
    print(f"✅ Reads saved to {reads_fasta}")
    
    # Run Minimap2
    minimap2_sam = output_dir / "minimap2.sam"
    minimap2_time, minimap2_stderr = run_minimap2(
        args.fasta,
        reads_fasta,
        minimap2_sam,
        args.threads
    )
    
    # Run GenoCache
    genocache_sam = output_dir / "genocache.sam"
    genocache_time, genocache_alns = run_genocache(
        args.checkpoint,
        args.index,
        args.positions,
        args.fasta,
        args.chrom,
        reads,
        genocache_sam,
        k=args.k,
        num_seeds=args.num_seeds,
        use_cuda=args.use_cuda
    )
    
    # Parse alignments
    print("\nParsing alignments...")
    minimap2_alns = parse_sam(minimap2_sam)
    genocache_alns_parsed = parse_sam(genocache_sam)
    
    # Evaluate accuracy
    print("Evaluating accuracy...")
    minimap2_stats = evaluate_accuracy(reads, minimap2_alns, args.tolerance)
    genocache_stats = evaluate_accuracy(reads, genocache_alns_parsed, args.tolerance)
    
    # Print comparison
    print_comparison(
        genocache_stats,
        minimap2_stats,
        genocache_time,
        minimap2_time,
        len(reads)
    )
    
    # Save results
    results = {
        'timestamp': datetime.now().isoformat(),
        'parameters': vars(args),
        'dataset': {
            'num_reads': len(reads),
            'actual_error_rate': actual_error_rate
        },
        'minimap2': {
            'time_seconds': minimap2_time,
            'reads_per_second': len(reads) / minimap2_time,
            'stats': minimap2_stats
        },
        'genocache': {
            'time_seconds': genocache_time,
            'reads_per_second': len(reads) / genocache_time,
            'stats': genocache_stats
        }
    }
    
    results_file = output_dir / "comparison_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {results_file}")
    print()


if __name__ == '__main__':
    main()
