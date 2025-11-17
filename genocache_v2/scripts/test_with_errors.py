#!/usr/bin/env python3
"""Test alignment with sequencing errors"""
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
    """Add random substitution errors"""
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

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint",required=True)
    p.add_argument("--index",required=True)
    p.add_argument("--positions",required=True)
    p.add_argument("--fasta",required=True)
    p.add_argument("--chrom",required=True)
    p.add_argument("--num-reads",type=int,default=1000)
    p.add_argument("--read-len",type=int,default=2000)
    p.add_argument("--k",type=int,default=10)
    p.add_argument("--error-rate",type=float,default=0.01)
    p.add_argument("--seed",type=int,default=42)
    args = p.parse_args()
    
    random.seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    print(f"Testing with {args.error_rate*100}% error rate")
    
    # Load everything
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    idx = faiss.read_index(args.index)
    positions = np.load(args.positions)
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    
    # Generate reads WITH ERRORS
    max_start = len(ref_seq) - args.read_len - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    reads = []
    for start in true_starts:
        clean_read = ref_seq[start:start + args.read_len]
        if 'N' not in clean_read:
            noisy_read = add_sequencing_errors(clean_read, args.error_rate)
            reads.append((start, noisy_read))
    
    print(f"Generated {len(reads)} reads with errors")
    
    # Align
    stats = {"total": 0, "mapped": 0, "no_candidate": 0, "correct_position": 0}
    tolerance = 1024  # Position tolerance
    
    for read_id, (true_pos, read_seq) in enumerate(tqdm(reads, desc="reads")):
        stats["total"] += 1
        candidates = []
        
        # 5 seeds per read
        seed_len = 512
        for i in range(5):
            offset = i * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            
            try:
                t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
                with torch.no_grad():
                    emb = model(t).cpu().detach().numpy()
                emb = np.ascontiguousarray(emb, dtype=np.float32)
                
                _, I = idx.search(emb, args.k)
                cands = positions[I[0]]
                for c in cands:
                    candidates.append(c - offset)
            except:
                continue
        
        if len(candidates) == 0:
            stats["no_candidate"] += 1
            continue
        
        # Pick most frequent
        unique, counts = np.unique(candidates, return_counts=True)
        best_candidate = unique[np.argmax(counts)]
        stats["mapped"] += 1
        
        # Check if correct
        if abs(best_candidate - true_pos) <= tolerance:
            stats["correct_position"] += 1
    
    # Results
    print(json.dumps({"stats": stats}, indent=2))
    print(f"\nAccuracy: {stats['correct_position']}/{stats['total']} = {100*stats['correct_position']/stats['total']:.1f}%")
    print(f"Mapping rate: {stats['mapped']}/{stats['total']} = {100*stats['mapped']/stats['total']:.1f}%")

if __name__ == "__main__":
    main()
