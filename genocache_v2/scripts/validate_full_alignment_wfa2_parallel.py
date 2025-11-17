#!/usr/bin/env python3
"""
FAISS -> WFA full-alignment validator
Writes SAM (and BAM if pysam installed) and JSON summary.
"""
import os, sys, argparse, time, tempfile, shlex, subprocess, json, re
from pathlib import Path
import numpy as np
from tqdm import tqdm

# local imports (run with PYTHONPATH=./scripts:.)
from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode

# ---------- WFA wrapper ----------
WFA_BIN = os.environ.get("WFA_BIN", "/home/nebius/tools/wfa2-lib/build/align_benchmark")

def call_wfa_align_fasta(query_seq, ref_seq_window, timeout=30):
    """Call align_benchmark with a two-sequence FASTA file and return parsed results."""
    if not os.path.exists(WFA_BIN):
        raise FileNotFoundError(f"WFA binary not found at {WFA_BIN}")
    tf = tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa")
    tf.write(">q\n" + query_seq + "\n")
    tf.write(">r\n" + ref_seq_window + "\n")
    tf.flush(); tf.close()
    cmd = f"{shlex.quote(WFA_BIN)} -a gap-affine-wfa -i {shlex.quote(tf.name)} --output-full /dev/stdout --quiet"
    try:
        p = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
        out = (p.stdout or "") + "\n" + (p.stderr or "")
    except subprocess.TimeoutExpired:
        out = ""
    finally:
        try: os.unlink(tf.name)
        except: pass
    score, ref_begin, ref_end, cigar = _parse_wfa_output(out)
    return score, ref_begin, ref_end, cigar, out

def _parse_wfa_output(out):
    cigar = None; score = None; ref_begin = None; ref_end = None
    m = re.search(r'(?i)refBegin[:=]?\s*(-?\d+)', out)
    if m: ref_begin = int(m.group(1))
    m = re.search(r'(?i)refEnd[:=]?\s*(-?\d+)', out)
    if m: ref_end = int(m.group(1))
    for line in out.splitlines():
        line = line.strip()
        m = re.match(r'^(-?\d+)\s*[\t ]\s*([0-9MIDNSHP=X]+)$', line)
        if m:
            score = int(m.group(1))
            cigar = m.group(2)
            break
    if cigar is None:
        m = re.search(r'(?i)(?:CIGAR|cigar)[:=]?\s*([0-9MIDNSHP=X]+)', out)
        if m: cigar = m.group(1)
    if score is None:
        m = re.search(r'(?i)score[:=]?\s*(-?\d+)', out)
        if m: score = int(m.group(1))
    if (score is None or ref_begin is None):
        nums = [int(x) for x in re.findall(r'(-?\d+)', out)]
        if len(nums) >= 1 and score is None:
            score = nums[0]
        if len(nums) >= 3 and (ref_begin is None or ref_end is None):
            if ref_begin is None: ref_begin = nums[1]
            if ref_end is None: ref_end = nums[2]
    return score, ref_begin, ref_end, cigar

# ---------- main pipeline ----------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--checkpoint', required=True)
    p.add_argument('--index', required=True)
    p.add_argument('--positions', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--num-reads', type=int, default=1000)
    p.add_argument('--read-len', type=int, default=2000)
    p.add_argument('--k', type=int, default=5)
    p.add_argument('--window', type=int, default=8192)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--out-bam', default="validation_wfa2.bam")
    p.add_argument('--use-cuda', action='store_true')
    return p.parse_args()

def write_sam_header(fout, chrom, chrom_len):
    fout.write(f"@HD\tVN:1.6\tSO:unsorted\n")
    fout.write(f"@SQ\tSN:{chrom}\tLN:{chrom_len}\n")
    fout.write("@PG\tID:genocache_validate\tPN:genocache_validate\tVN:0.1\n")

def main():
    args = parse_args()
    t_start = time.time()

    device = "cuda" if (args.use_cuda and __import__("torch").cuda.is_available()) else "cpu"
    model, meta = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    seed_len = model.input_len

    ref = ReferenceGenome(args.fasta)
    seq = ref.sequences[args.chrom]
    chrom_len = len(seq)
    positions = np.load(args.positions)
    import faiss
    idx = faiss.read_index(args.index)
    print("Index loaded:", idx.ntotal, "vectors")

    rng = np.random.default_rng(args.seed)
    starts = []
    while len(starts) < args.num_reads:
        p = int(rng.integers(0, max(1, chrom_len - args.read_len)))
        r = seq[p:p+args.read_len]
        if r.count('N') > args.read_len * 0.1:
            continue
        starts.append(p)

    out_sam = args.out_bam.replace(".bam", ".sam")
    fout = open(out_sam, "w")
    write_sam_header(fout, args.chrom, chrom_len)

    stats = {'total':0,'mapped':0,'unmapped':0,'no_candidate':0,'wfa_fail':0}
    per_read = []

    for i, start in enumerate(tqdm(starts, desc="reads")):
        stats['total'] += 1
        read_seq = seq[start:start+args.read_len]
        qname = f"read_{i}_pos{start}"
        best_hit = {'score': None, 'ref_start': None, 'ref_end': None, 'cigar': None, 'raw': None, 'candidate_pos': None}
        seed_offsets = [0, 128, 256, 384, 512]
        found_candidates = False

        for offset in seed_offsets:
            if offset + seed_len > len(read_seq):
                continue
            seed_seq = read_seq[offset:offset+seed_len]
            import torch
            seed_tensor = one_hot_encode(seed_seq)
            import numpy as _np
            import torch as _torch
            # accept numpy array, torch tensor, or Python sequence
            if isinstance(seed_tensor, _np.ndarray):
                seed_tensor = _torch.from_numpy(seed_tensor)
            elif hasattr(seed_tensor, 'detach') and hasattr(seed_tensor, 'unsqueeze'):
                # already a torch tensor
                pass
            else:
                seed_tensor = _torch.tensor(seed_tensor)
            t = seed_tensor.unsqueeze(0).to(device).float()
            # Embed seed safely
            emb_t = model(t).detach().cpu()
            try:
                emb = emb_t.numpy()
            except Exception:
                emb = _np.array(emb_t)
            D, I = idx.search(emb, args.k)
            cand_idx = I[0]
            if cand_idx.size == 0:
                continue
            found_candidates = True

            for c in cand_idx:
                if c < len(positions):
                    cand_pos = int(positions[c])
                else:
                    cand_pos = int(positions[c % len(positions)])
                w = args.window
                ref_s = max(0, cand_pos - w//2)
                ref_e = min(chrom_len, ref_s + w)
                ref_window_seq = seq[ref_s:ref_e]
                score, rb, re, cigar, raw = call_wfa_align_fasta(seed_seq, ref_window_seq, timeout=20)
                if rb is not None:
                    genome_rb = ref_s + rb
                    genome_re = ref_s + (re if re is not None else (rb + len(seed_seq)))
                else:
                    genome_rb = cand_pos
                    genome_re = cand_pos + len(seed_seq)
                if score is None:
                    stats['wfa_fail'] += 1
                    continue
                if (best_hit['score'] is None) or (score > best_hit['score']):
                    best_hit.update({
                        'score': score, 'ref_start': genome_rb, 'ref_end': genome_re,
                        'cigar': cigar or f"{len(seed_seq)}M", 'raw': raw, 'candidate_pos': cand_pos,
                        'seed_seq': seed_seq
                    })
            if best_hit['score'] is not None and best_hit['score'] > -9999:
                break

        if not found_candidates:
            stats['no_candidate'] += 1
            stats['unmapped'] += 1
            fout.write(f"{qname}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
            per_read.append({'qname':qname,'mapped':False})
            continue

        if best_hit['score'] is None:
            stats['unmapped'] += 1
            fout.write(f"{qname}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
            per_read.append({'qname':qname,'mapped':False})
            continue

        pos1 = max(1, int(best_hit['ref_start']) + 1)
        mapq = 255
        cigar = best_hit['cigar'] or f"{args.read_len}M"
        seq_to_write = best_hit.get('seed_seq', read_seq)
        fout.write(f"{qname}\t0\t{args.chrom}\t{pos1}\t{mapq}\t{cigar}\t*\t0\t0\t{seq_to_write}\t*\n")
        stats['mapped'] += 1
        per_read.append({'qname':qname,'mapped':True,'pos':int(best_hit['ref_start']),'score':best_hit['score'],'cigar':cigar})

    fout.close()

    out_bam = args.out_bam
    try:
        import pysam
        pysam.view("-bS", out_sam, "-o", out_bam, catch_stdout=False)
        pysam.index(out_bam)
    except Exception:
        out_bam = None

    summary = {
        'args': vars(args),
        'time_s': time.time() - t_start,
        'stats': stats,
        'num_reads': len(starts),
        'mapped': stats['mapped'],
    }
    with open("validation_wfa2_summary.json", "w") as jf:
        json.dump(summary, jf, indent=2)

    print("Done. Summary saved to validation_wfa2_summary.json")
    print(summary)

if __name__ == "__main__":
    main()
