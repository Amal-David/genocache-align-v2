#!/usr/bin/env python3
"""
align_full_reads_with_outputs.py

Patched alignment runner that:
- loads trained ImprovedCNN checkpoint
- loads FAISS index and reference positions
- generates synthetic reads (or loads from fasta if implemented)
- finds candidate reference positions via FAISS
- (optionally) calls an external WFA binary to verify and produce CIGARs
- writes: per-read CSV, sorted BAM (pysam), summary JSON

Usage (example):
python scripts/align_full_reads_with_outputs.py \
  --checkpoint /home/nebius/work/genocache_checkpoints/improved_cnn_best.pt \
  --index indexes/index_ivfpq.faiss \
  --positions data/reference_encodings/ref_positions_20251111_163637.npy \
  --fasta references/GCF_000001405.40_GRCh38.p14_chr22.fna \
  --chrom NC_000022.11 \
  --num-reads 1000 \
  --read-len 2000 \
  --k 10 \
  --window 8192 \
  --out-prefix align_chr22_1000

Notes:
- If WFA_BIN environment variable is set and points to a compatible align_benchmark, this script will attempt to call it to compute exact CIGAR strings. If not present, it will still write BAM + CSV with placeholder CIGARs and mark wfa_status accordingly.
- Requires: numpy, faiss, torch, pysam, tqdm
"""

import argparse
import json
import os
import sys
import time
import tempfile
import subprocess
from pathlib import Path
from datetime import datetime

import numpy as np
import faiss
import torch

# local imports: rely on PYTHONPATH=./scripts:.
from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode

try:
    import pysam
except Exception as e:
    pysam = None


def load_index(idx_path):
    idx = faiss.read_index(str(idx_path))
    return idx


def load_positions(pos_path):
    return np.load(pos_path)


def generate_reads_from_reference(fasta_path, chrom, num_reads, read_len, seed=None):
    ref = ReferenceGenome(fasta_path)
    seq = ref.sequences[chrom]
    seq_len = len(seq)
    rng = np.random.RandomState(seed)
    reads = []
    positions = []
    while len(reads) < num_reads:
        pos = int(rng.randint(0, seq_len - read_len))
        s = seq[pos:pos + read_len]
        if s.count('N') > int(read_len * 0.1):
            continue
        reads.append(s)
        positions.append(pos)
    return reads, positions


def maybe_run_wfa(wfa_bin, query_seq, ref_seq, algo='gap-affine-wfa', penalties=None, timeout=5.0):
    """
    Attempt to run align_benchmark (wfa) on a single pair. Because align_benchmark requires an input file, write a tiny FASTA and call it.
    Returns: (wfa_status, cigar, wfa_time)
    If wfa_bin is None or call fails, returns ('SKIPPED','*',0.0)
    """
    if not wfa_bin:
        return 'WFA_MISSING', '*', 0.0

    # write temporary FASTA with two sequences: pattern and text using simple format expected by align_benchmark
    with tempfile.TemporaryDirectory() as td:
        in_fasta = Path(td) / 'pairs.fa'
        out_file = Path(td) / 'wfa_out.txt'
        # align_benchmark expects an input; format varies, but many builds accept simple two-sequence FASTA pairs per record
        with open(in_fasta, 'w') as f:
            f.write(f">q\n{query_seq}\n>t\n{ref_seq}\n")
        cmd = [wfa_bin, '-a', 'gap-affine-wfa', '-i', str(in_fasta), '--output', str(out_file), '--quiet']
        if penalties:
            cmd += ['--affine-penalties', penalties]
        t0 = time.time()
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
            t1 = time.time()
            wtime = t1 - t0
            # parse out_file for CIGAR (best-effort)
            cigar = '*'
            if out_file.exists():
                text = out_file.read_text(errors='ignore')
                # Try to find a CIGAR-like string or alignment line; fallback to '*'
                # This is heuristic — adapt if your align_benchmark output format differs
                for line in text.splitlines():
                    if '\t' in line and ('M' in line or 'I' in line or 'D' in line or 'X' in line):
                        cigar = line.split('\t')[-1].strip()
                        break
            return 'OK', cigar, wtime
        except subprocess.CalledProcessError:
            return 'WFA_FAIL', '*', 0.0
        except Exception:
            return 'WFA_ERROR', '*', 0.0


def write_bam_and_index(out_bam_path, records, ref_name, ref_len):
    if pysam is None:
        raise RuntimeError('pysam is required to write BAM. Please install pysam in your environment.')
    # create header
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'SN': ref_name, 'LN': int(ref_len)}]}
    # write unsorted BAM to temp then sort
    tmp_bam = str(out_bam_path) + '.unsorted.bam'
    with pysam.AlignmentFile(tmp_bam, 'wb', header=header) as outf:
        for rec in records:
            outf.write(rec)
    sorted_bam = str(out_bam_path)
    pysam.sort('-o', sorted_bam, tmp_bam)
    pysam.index(sorted_bam)
    os.remove(tmp_bam)


def make_aligned_segment(read_name, read_seq, ref_name, ref_pos0, cigar, mapq=60, flag=0):
    # create pysam.AlignedSegment
    a = pysam.AlignedSegment()
    a.query_name = read_name
    a.flag = flag
    
    if flag & 4:  # Unmapped
        a.reference_id = -1
        a.reference_start = -1
        a.mapping_quality = 0
        a.cigarstring = None
    else:  # Mapped
        a.reference_id = 0
        a.reference_start = int(ref_pos0)
        a.mapping_quality = mapq
        
        # Use actual CIGAR from WFA
        if cigar and cigar != '*':
            a.cigarstring = cigar
        else:
            a.cigarstring = None
    
    a.query_sequence = read_seq
    # For simplicity, set mate fields to default
    return a


def main():
    parser = argparse.ArgumentParser(description='Align reads using GenoCache-Align pipeline and write outputs (BAM + CSV + JSON summary)')
    parser.add_argument('--checkpoint', required=True)
    parser.add_argument('--index', required=True)
    parser.add_argument('--positions', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--chrom', required=True)
    parser.add_argument('--num-reads', type=int, default=1000)
    parser.add_argument('--read-len', type=int, default=2000)
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--window', type=int, default=8192)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--out-prefix', type=str, required=True)
    parser.add_argument('--use-cuda', action='store_true')
    args = parser.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # outputs
    per_read_csv = out_prefix.parent / (out_prefix.name + '_per_read.csv')
    out_bam = out_prefix.parent / (out_prefix.name + '.bam')
    summary_json = out_prefix.parent / (out_prefix.name + '_summary.json')

    device = 'cuda' if args.use_cuda and torch.cuda.is_available() else 'cpu'

    print('Loading model...')
    model, metadata = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    seed_len = model.input_len

    print('Loading index...')
    idx = load_index(args.index)
    print(f'Index loaded: {idx.ntotal} vectors')

    print('Loading positions...')
    positions = load_positions(args.positions)

    # reference object for extracting ref substrings for WFA if needed
    ref = ReferenceGenome(args.fasta)
    ref_seq = ref.sequences[args.chrom]
    ref_len = len(ref_seq)

    # generate reads
    print('Generating synthetic reads...')
    reads, true_positions = generate_reads_from_reference(args.fasta, args.chrom, args.num_reads, args.read_len, seed=args.seed)

    # prepare outputs
    per_read_f = open(per_read_csv, 'w')
    per_read_f.write('read_id,true_pos,chosen_candidate,seed_offset,success,wfa_status,wfa_time_s,cigar\n')

    records = []
    stats = {'total': 0, 'mapped': 0, 'unmapped': 0, 'no_candidate': 0, 'wfa_fail': 0}

    # optional WFA binary
    wfa_bin = os.environ.get('WFA_BIN', None)
    if wfa_bin:
        print('WFA binary:', wfa_bin)
    else:
        print('WFA binary not set; WFA step will be skipped (wfa_status=SKIPPED)')

    t_start = time.time()

    for i, (read_seq, true_pos) in enumerate(zip(reads, true_positions)):
        stats['total'] += 1
        read_name = f'read_{i}_pos{true_pos}'
        # iterate seeds across read to find candidate
        found = False
        chosen = None
        chosen_cigar = '*'
        chosen_wfa_time = 0.0
        chosen_wfa_status = 'SKIPPED' if not wfa_bin else 'UNTRIED'
        seed_offsets = list(range(0, args.read_len - seed_len + 1, max(1, seed_len//4)))
        # limit number of seeds per read for speed
        seed_offsets = seed_offsets[0:5]

        for offset in seed_offsets:
            seed_seq = read_seq[offset:offset + seed_len]
            if len(seed_seq) != seed_len:
                continue
            t = torch.tensor(one_hot_encode(seed_seq)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(t).cpu().numpy()
            # search - ensure contiguous float32 array
            emb = np.ascontiguousarray(emb, dtype=np.float32)
            D, I = idx.search(emb, args.k)
            cand_idx = I[0]
            if cand_idx.size == 0:
                continue
            # map candidate indices to positions
            cand_positions = positions[cand_idx]
            # pick top candidate
            chosen_idx = cand_idx[0]
            chosen_pos = int(cand_positions[0])
            # optionally run WFA on a small window around candidate
            if wfa_bin:
                ref_start = max(0, chosen_pos - 100)
                ref_end = min(ref_len, chosen_pos + args.read_len + 100)
                ref_subseq = ref_seq[ref_start:ref_end]
                wfa_status, cigar, wtime = maybe_run_wfa(wfa_bin, read_seq, ref_subseq)
                chosen_wfa_status = wfa_status
                chosen_cigar = cigar
                chosen_wfa_time = wtime
                if wfa_status == 'OK':
                    chosen = chosen_pos if ref_start == 0 else chosen_pos
                    found = True
                    break
                else:
                    # accept candidate even if WFA failed (best-effort)
                    chosen = chosen_pos
                    found = True
                    if wfa_status.startswith('WFA_FAIL'):
                        stats['wfa_fail'] += 1
                    break
            else:
                chosen = chosen_pos
                found = True
                break

        if found and chosen is not None:
            stats['mapped'] += 1
            success = 1
            # create aligned record (pysam) if pysam available
            if pysam is not None:
                # simple CIGAR: assume match of read_len
                if chosen_cigar != '*':
                    # Valid alignment
                    print(f"Creating record: {read_name}, pos={chosen}, cigar={chosen_cigar}, read_len={len(read_seq)}")
                    rec = make_aligned_segment(read_name, read_seq, args.chrom, chosen, chosen_cigar)
                else:
                    # Unmapped read - set flag=4 and no position/cigar
                    rec = make_aligned_segment(read_name, read_seq, args.chrom, 0, None, mapq=0, flag=4)
                records.append(rec)
            per_read_f.write(f"{read_name},{true_pos},{chosen},{offset},{success},{chosen_wfa_status},{chosen_wfa_time},{chosen_cigar}\n")
        else:
            stats['unmapped'] += 1
            per_read_f.write(f"{read_name},{true_pos},-1,-1,0,NO_CAND,0,*\n")

    t_end = time.time()
    per_read_f.close()

    # write BAM if pysam available
    if pysam is not None and records:
        print('Writing BAM...')
        write_bam_and_index(out_bam, records, args.chrom, ref_len)
        print('BAM written:', out_bam)
    else:
        print('pysam not available or no records; skipping BAM write')

    summary = {
        'args': vars(args),
        'time_s': t_end - t_start,
        'stats': stats,
        'num_reads': stats['total'],
        'mapped': stats['mapped']
    }
    with open(summary_json, 'w') as f:
        json.dump(summary, f, indent=2)

    print('Done. Summary saved to', summary_json)


if __name__ == '__main__':
    main()
