#!/usr/bin/env python3
import argparse, json, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
import sys
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode

def load_fasta(fasta_path):
    """Load FASTA file into dictionary"""
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
    p.add_argument("--num-reads",type=int,default=100)
    p.add_argument("--read-len",type=int,default=2000)
    p.add_argument("--k",type=int,default=5)
    p.add_argument("--seed",type=int,default=42)
    args = p.parse_args()
    
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    # Load model
    print(f"Loading model...")
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    
    # Load FAISS
    print(f"Loading FAISS index...")
    idx = faiss.read_index(args.index)
    positions = np.load(args.positions)
    print(f"Index: {idx.ntotal:,} vectors")
    
    # Load reference
    print(f"Loading reference...")
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    print(f"Reference: {len(ref_seq):,} bp")
    
    # Generate test reads
    print(f"Generating {args.num_reads} reads...")
    max_start = len(ref_seq) - args.read_len - 10
    true_starts = np.random.randint(10, max_start, args.num_reads)
    reads = []
    for start in true_starts:
        read_seq = ref_seq[start:start + args.read_len]
        if 'N' not in read_seq:
            reads.append((start, read_seq))
    print(f"Valid reads: {len(reads)}")
    
    # Align
    stats = {"total": 0, "mapped": 0, "no_candidate": 0}
    
    print(f"Aligning...")
    for read_id, (true_pos, read_seq) in enumerate(tqdm(reads, desc="reads")):
        stats["total"] += 1
        candidates = []
        
        # 5 seeds per read
        seed_len = 512
        for i in range(5):
            offset = i * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            
            # Encode
            t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
            with torch.no_grad():
                emb = model(t).cpu().detach().numpy()
            emb = np.ascontiguousarray(emb, dtype=np.float32)
            
            # Search
            _, I = idx.search(emb, args.k)
            cands = positions[I[0]]
            for c in cands:
                candidates.append(c - offset)
        
        if len(candidates) == 0:
            stats["no_candidate"] += 1
            continue
        
        # Pick most frequent
        unique, counts = np.unique(candidates, return_counts=True)
        stats["mapped"] += 1
    
    # Results
    print(json.dumps({"stats": stats}, indent=2))

if __name__ == "__main__":
    main()
