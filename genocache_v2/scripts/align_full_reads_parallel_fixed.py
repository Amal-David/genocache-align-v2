#!/usr/bin/env python3
import argparse, time, json, subprocess, os, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode, ReferenceGenome

def main(args):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    
    idx = faiss.read_index(args.index)
    positions = np.load(args.positions)
    ref = ReferenceGenome(args.fasta)
    ref_seq = ref.get_sequence(args.chrom)
    
    # Generate reads
    max_start = len(ref_seq) - args.read_len - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    reads = [(s, ref_seq[s:s+args.read_len]) for s in true_starts if 'N' not in ref_seq[s:s+args.read_len]]
    
    stats = {"total": 0, "mapped": 0, "unmapped": 0, "no_candidate": 0, "wfa_fail": 0}
    
    for read_id, (true_pos, read_seq) in enumerate(tqdm(reads, desc="reads")):
        stats["total"] += 1
        candidates = []
        
        # Extract seeds
        seed_len = 512
        num_seeds = min(5, (len(read_seq) - seed_len) // 100 + 1)
        for i in range(num_seeds):
            offset = i * ((len(read_seq) - seed_len) // max(1, num_seeds - 1))
            seed = read_seq[offset:offset+seed_len]
            
            # Encode with proper numpy conversion
            t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(t).cpu().detach().numpy()
            
            # CRITICAL: Ensure proper numpy format for FAISS
            emb = np.ascontiguousarray(emb, dtype=np.float32)
            
            _, I = idx.search(emb, args.k)
            cands = positions[I[0]]
            if len(cands)==0:
                continue
            for c in cands:
                candidates.append(c - offset)
        
        if len(candidates)==0:
            stats["no_candidate"]+=1
            continue
        
        # Pick most frequent
        unique, counts = np.unique(candidates, return_counts=True)
        best_candidate = unique[np.argmax(counts)]
        stats["mapped"] += 1
    
    print(json.dumps({"stats": stats}, indent=2))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",required=True)
    p.add_argument("--index",required=True)
    p.add_argument("--positions",required=True)
    p.add_argument("--fasta",required=True)
    p.add_argument("--chrom",required=True)
    p.add_argument("--num-reads",type=int,default=100)
    p.add_argument("--read-len",type=int,default=2000)
    p.add_argument("--k",type=int,default=5)
    p.add_argument("--window",type=int,default=8192)
    p.add_argument("--seed",type=int,default=42)
    p.add_argument("--wfa-bin",default="WFA2-lib/build/align_benchmark")
    args = p.parse_args()
    np.random.seed(args.seed)
    main(args)
