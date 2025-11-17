#!/usr/bin/env python3
"""Test with errors + clustered voting"""
import random
import argparse, json, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
import sys
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode

def add_sequencing_errors(seq, error_rate=0.01):
    seq = list(seq)
    bases = ['A', 'C', 'G', 'T']
    for i in range(len(seq)):
        if random.random() < error_rate:
            seq[i] = random.choice([b for b in bases if b != seq[i]])
    return ''.join(seq)

def load_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences

def cluster_positions(positions, window=100):
    if len(positions) == 0:
        return []
    positions = sorted(positions)
    clusters = []
    current_cluster = [positions[0]]
    for pos in positions[1:]:
        if pos - current_cluster[-1] <= window:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)
    return clusters

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",required=True)
    p.add_argument("--index",required=True)
    p.add_argument("--positions",required=True)
    p.add_argument("--fasta",required=True)
    p.add_argument("--chrom",required=True)
    p.add_argument("--num-reads",type=int,default=1000)
    p.add_argument("--error-rate",type=float,default=0.01)
    p.add_argument("--cluster-window",type=int,default=100)
    p.add_argument("--tolerance",type=int,default=1024)
    args = p.parse_args()
    
    random.seed(42)
    np.random.seed(42)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    print(f"Testing with {args.error_rate*100}% error rate + clustered voting")
    
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    index = faiss.read_index(args.index)
    ref_positions = np.load(args.positions)
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    
    max_start = len(ref_seq) - 2000 - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    reads = []
    for start in true_starts:
        clean_read = ref_seq[start:start + 2000]
        if 'N' not in clean_read:
            noisy_read = add_sequencing_errors(clean_read, args.error_rate)
            reads.append((start, noisy_read))
    
    print(f"Generated {len(reads)} reads with {args.error_rate*100}% errors")
    
    stats = {"total": 0, "mapped": 0, "correct": 0}
    
    for read_id, (true_pos, read_seq) in enumerate(tqdm(reads, desc="reads")):
        stats["total"] += 1
        candidates = []
        
        seed_len = 512
        for i in range(5):
            offset = i * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            
            try:
                t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
                with torch.no_grad():
                    emb = model(t).cpu().detach().numpy()
                emb = np.ascontiguousarray(emb, dtype=np.float32)
                
                D, I = index.search(emb, k=10)
                for c in ref_positions[I[0]]:
                    candidates.append(c - offset)
            except:
                continue
        
        if not candidates:
            continue
        
        clusters = cluster_positions(candidates, window=args.cluster_window)
        best_cluster = max(clusters, key=len)
        best_position = int(np.median(best_cluster))
        
        stats["mapped"] += 1
        if abs(best_position - true_pos) <= args.tolerance:
            stats["correct"] += 1
    
    print(json.dumps({"stats": stats}, indent=2))
    print(f"\nAccuracy: {stats['correct']}/{stats['total']} = {100*stats['correct']/stats['total']:.1f}%")
    print(f"Mapping rate: {stats['mapped']}/{stats['total']} = {100*stats['mapped']/stats['total']:.1f}%")

if __name__ == "__main__":
    main()
