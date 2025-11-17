#!/usr/bin/env python3
"""
End-to-end alignment pipeline using neural encoder + FAISS index
"""
import os, sys, time
import numpy as np
import torch
import torch.nn as nn
import faiss
from Bio import SeqIO
from collections import defaultdict

# Config
ENCODER_PATH = "nal_encoder_epoch1.pt"
FAISS_INDEX = "faiss_index_cpu.ivf"
REF_POSITIONS = "ref_positions.npy"
REF_FA = "chr22_1kb.fa"
QUERY_FA = "genocache_data/reads_chr22_synth_1kb_500.fa"
TRUTH_TSV = "genocache_data/reads_chr22_truth.tsv"

SEED_LEN = 256
TOP_K = 10
TOLERANCE = 64  # bases

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

class TinyNAL(nn.Module):
    def __init__(self, out_dim=128):
        super().__init__()
        self.in_proj = nn.Conv1d(4, 64, kernel_size=7, padding=3)
        self.conv = nn.Sequential(
            nn.ReLU(),
            nn.Conv1d(64, 128, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.proj = nn.Linear(128, out_dim)
    
    def forward(self, x):
        x = x.permute(0,2,1)
        x = self.in_proj(x)
        x = self.conv(x)
        x = x.view(x.shape[0], -1)
        x = self.proj(x)
        return nn.functional.normalize(x, dim=-1)

def seq_to_onehot(seq):
    import random
    arr = np.zeros((len(seq),4), dtype=np.float32)
    for i,ch in enumerate(seq):
        if ch.upper() == "A": arr[i,0]=1
        elif ch.upper()=="C": arr[i,1]=1
        elif ch.upper()=="G": arr[i,2]=1
        elif ch.upper()=="T": arr[i,3]=1
        elif ch.upper()=="N": arr[i, random.randrange(4)] = 0.25  # uniform
        else: arr[i, random.randrange(4)] = 1.0
    return arr

def load_truth_positions(truth_file):
    truth = {}
    with open(truth_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                read_id, start, end = parts[0], int(parts[1]), int(parts[2])
                truth[read_id] = start
    return truth

def main():
    print("=" * 60)
    print("GenoCache Alignment Pipeline")
    print("=" * 60)
    
    # Load encoder
    print(f"\n[1/6] Loading encoder: {ENCODER_PATH}")
    model = TinyNAL().to(DEVICE)
    model.load_state_dict(torch.load(ENCODER_PATH, map_location=DEVICE))
    model.eval()
    print(f"  ✓ Encoder loaded on {DEVICE}")
    
    # Load FAISS index
    print(f"\n[2/6] Loading FAISS index: {FAISS_INDEX}")
    index = faiss.read_index(FAISS_INDEX)
    if hasattr(index, 'nprobe'):
        index.nprobe = 16  # increase search thoroughness
    print(f"  ✓ Index loaded: {index.ntotal} vectors")
    
    # Load reference positions
    print(f"\n[3/6] Loading reference positions")
    ref_pos = np.load(REF_POSITIONS)
    print(f"  ✓ Loaded {len(ref_pos)} position mappings")
    
    # Load truth
    print(f"\n[4/6] Loading ground truth: {TRUTH_TSV}")
    truth = load_truth_positions(TRUTH_TSV)
    print(f"  ✓ Loaded truth for {len(truth)} reads")
    
    # Load and encode query reads
    print(f"\n[5/6] Loading query reads: {QUERY_FA}")
    queries = []
    query_ids = []
    for record in SeqIO.parse(QUERY_FA, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < SEED_LEN:
            continue
        # Take first SEED_LEN bases
        seq = seq[:SEED_LEN]
        queries.append(seq)
        query_ids.append(record.id)
    
    print(f"  ✓ Loaded {len(queries)} queries")
    print(f"\n  Encoding queries...")
    
    batch_size = 128
    all_vecs = []
    with torch.no_grad():
        for i in range(0, len(queries), batch_size):
            batch_seqs = queries[i:i+batch_size]
            batch_arr = np.stack([seq_to_onehot(s) for s in batch_seqs])
            batch_tensor = torch.from_numpy(batch_arr).to(DEVICE)
            vecs = model(batch_tensor).cpu().numpy()
            all_vecs.append(vecs)
    
    query_vecs = np.vstack(all_vecs).astype('float32')
    print(f"  ✓ Encoded {len(query_vecs)} query vectors")
    
    # Search index
    print(f"\n[6/6] Searching index (top-{TOP_K})...")
    start_time = time.time()
    D, I = index.search(query_vecs, TOP_K)
    search_time = time.time() - start_time
    print(f"  ✓ Search completed in {search_time:.2f}s ({len(queries)/search_time:.1f} queries/sec)")
    
    # Evaluate results
    print(f"\n{'='*60}")
    print("Results")
    print("=" * 60)
    
    correct_by_k = defaultdict(int)
    total = len(query_ids)
    
    for qi, qid in enumerate(query_ids):
        if qid not in truth:
            continue
        
        true_pos = truth[qid]
        top_indices = I[qi]
        top_positions = ref_pos[top_indices]
        
        for k in [1, 5, 10]:
            if k > TOP_K:
                continue
            candidates = top_positions[:k]
            if np.any(np.abs(candidates - true_pos) <= TOLERANCE):
                correct_by_k[k] += 1
    
    print(f"\nRecall @ tolerance={TOLERANCE}bp:")
    for k in sorted(correct_by_k.keys()):
        recall = correct_by_k[k] / total
        print(f"  Recall@{k:2d}: {correct_by_k[k]:4d}/{total} = {recall*100:5.2f}%")
    
    # Show examples
    print(f"\nExample alignments (first 5):")
    for qi in range(min(5, len(query_ids))):
        qid = query_ids[qi]
        if qid not in truth:
            continue
        true_pos = truth[qid]
        pred_pos = ref_pos[I[qi, 0]]
        error = abs(pred_pos - true_pos)
        status = "✓" if error <= TOLERANCE else "✗"
        print(f"  {status} {qid}: true={true_pos:9d}  pred={pred_pos:9d}  error={error:5d}bp")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
