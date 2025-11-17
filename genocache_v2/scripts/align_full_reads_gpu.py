#!/usr/bin/env python3
"""GPU-Optimized Full Read Aligner with FAISS on GPU"""

import argparse, time, json, subprocess, os, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode, ReferenceGenome

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--num-reads", type=int, default=100)
    p.add_argument("--read-len", type=int, default=2000)
    p.add_argument("--k", type=int, default=5)
    p.add_argument("--window", type=int, default=8192)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"🚀 Device: {device}")
    
    # Load model
    model, _ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    
    # Load FAISS and move to GPU
    print(f"📂 Loading FAISS...")
    idx_cpu = faiss.read_index(args.index)
    print(f"   Vectors: {idx_cpu.ntotal:,}")
    
    if device == "cuda":
        print(f"⚡ Moving FAISS to GPU...")
        res = faiss.StandardGpuResources()
        idx = faiss.index_cpu_to_gpu(res, 0, idx_cpu)
        print(f"   ✅ FAISS on GPU")
    else:
        idx = idx_cpu
    
    positions = np.load(args.positions)
    ref = ReferenceGenome(args.fasta)
    ref_seq = ref.get_sequence(args.chrom)
    print(f"🧬 Ref: {len(ref_seq):,} bp")
    
    # Generate reads
    max_start = len(ref_seq) - args.read_len - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    reads = [(s, ref_seq[s:s+args.read_len]) for s in true_starts if 'N' not in ref_seq[s:s+args.read_len]]
    print(f"🎲 Reads: {len(reads)}")
    
    stats = {"total": 0, "mapped": 0, "unmapped": 0, "no_candidate": 0}
    results = []
    
    print(f"🔍 Aligning...")
    start = time.time()
    
    for rid, (true_pos, read_seq) in enumerate(tqdm(reads, desc="reads")):
        stats["total"] += 1
        candidates = []
        
        # 5 seeds per read
        for i in range(5):
            offset = i * ((len(read_seq) - 512) // 4)
            seed = read_seq[offset:offset+512]
            t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(t).cpu().numpy()
            D, I = idx.search(emb, args.k)
            candidates.extend(positions[I[0]] - offset)
        
        if not candidates:
            stats["no_candidate"] += 1
            continue
        
        # Pick most frequent
        unique, counts = np.unique(candidates, return_counts=True)
        best = unique[np.argmax(counts)]
        stats["mapped"] += 1
        results.append({"read_id": rid, "true_pos": int(true_pos), "best": int(best)})
    
    elapsed = time.time() - start
    
    print(f"\n{'='*60}")
    print(f"⏱️  Time: {elapsed:.2f}s")
    print(f"🚀 Throughput: {stats['total']/elapsed:.2f} reads/sec")
    print(json.dumps({"stats": stats}, indent=2))
    print(f"{'='*60}")
    
    # Save
    with open(f"align_{args.num_reads}reads_gpu_summary.json", 'w') as f:
        json.dump({"stats": stats, "time_s": elapsed, "throughput": stats['total']/elapsed}, f, indent=2)
    
    with open(f"align_{args.num_reads}reads_gpu_per_read.csv", 'w') as f:
        f.write("read_id,true_pos,best_candidate\n")
        for r in results:
            f.write(f"{r['read_id']},{r['true_pos']},{r['best']}\n")

if __name__ == "__main__":
    main()
