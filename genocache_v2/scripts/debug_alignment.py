#!/usr/bin/env python3
"""Debug alignment to see what's happening"""
import argparse, json, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
import sys
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode

def load_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",required=True)
    p.add_argument("--index",required=True)
    p.add_argument("--positions",required=True)
    p.add_argument("--fasta",required=True)
    p.add_argument("--chrom",required=True)
    p.add_argument("--num-reads",type=int,default=10)
    args = p.parse_args()
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    idx = faiss.read_index(args.index)
    positions = np.load(args.positions)
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    
    print(f"Index size: {idx.ntotal:,} vectors")
    print(f"Reference length: {len(ref_seq):,} bp")
    print(f"Position array size: {len(positions):,}")
    
    # Generate perfect reads
    np.random.seed(42)
    max_start = len(ref_seq) - 2000 - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    
    for i, true_pos in enumerate(true_starts[:args.num_reads]):
        read_seq = ref_seq[true_pos:true_pos + 2000]
        if 'N' in read_seq:
            continue
        
        print(f"\n{'='*60}")
        print(f"Read {i}: True position = {true_pos}")
        
        candidates = []
        seed_len = 512
        
        for seed_idx in range(5):
            offset = seed_idx * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            
            t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(t).cpu().detach().numpy()
            emb = np.ascontiguousarray(emb, dtype=np.float32)
            
            D, I = idx.search(emb, k=10)
            
            print(f"  Seed {seed_idx} (offset={offset}):")
            print(f"    Top-5 scores: {D[0][:5]}")
            print(f"    Top-5 positions: {positions[I[0]][:5]}")
            print(f"    Adjusted positions: {positions[I[0]][:5] - offset}")
            print(f"    Distance to true: {np.abs(positions[I[0]][:5] - offset - true_pos)}")
            
            for c in positions[I[0]]:
                candidates.append(c - offset)
        
        # Voting
        unique, counts = np.unique(candidates, return_counts=True)
        best_candidate = unique[np.argmax(counts)]
        
        print(f"\n  All candidates: {len(candidates)}")
        print(f"  Unique positions: {len(unique)}")
        print(f"  Top 5 voted positions:")
        top_indices = np.argsort(counts)[-5:][::-1]
        for idx in top_indices:
            print(f"    Position {unique[idx]}: {counts[idx]} votes, error = {abs(unique[idx] - true_pos)} bp")
        print(f"  Best candidate: {best_candidate} (votes: {np.max(counts)})")
        print(f"  True position: {true_pos}")
        print(f"  Error: {abs(best_candidate - true_pos)} bp")
        print(f"  Correct: {abs(best_candidate - true_pos) <= 1024}")

if __name__ == "__main__":
    main()
