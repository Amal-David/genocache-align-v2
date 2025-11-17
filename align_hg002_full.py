#!/usr/bin/env python3
"""
Align HG002 reads to full GRCh38 genome
"""
import os, sys, time, gzip
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import faiss
from collections import defaultdict
from Bio import SeqIO

# Config
ENCODER_PATH = "nal_encoder_best.pt"
FAISS_INDEX = "faiss_grch38.idx"
REF_POSITIONS = "grch38_positions.npy"
REF_CHROMOSOMES = "grch38_chromosomes.npy"

SEED_LEN = 256
EMB_DIM = 256
TOP_K = 20
BATCH_SIZE = 512

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

class ImprovedNAL(nn.Module):
    def __init__(self, out_dim=EMB_DIM):
        super().__init__()
        self.in_proj = nn.Conv1d(4, 128, kernel_size=9, padding=4)
        self.bn1 = nn.BatchNorm1d(128)
        self.conv1 = nn.Conv1d(128, 256, kernel_size=7, padding=3)
        self.bn2 = nn.BatchNorm1d(256)
        self.conv2 = nn.Conv1d(256, 256, kernel_size=7, padding=3)
        self.bn3 = nn.BatchNorm1d(256)
        self.conv3 = nn.Conv1d(256, 512, kernel_size=5, padding=2)
        self.bn4 = nn.BatchNorm1d(512)
        self.conv4 = nn.Conv1d(512, 512, kernel_size=5, padding=2)
        self.bn5 = nn.BatchNorm1d(512)
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.proj = nn.Sequential(nn.Linear(512, out_dim), nn.LayerNorm(out_dim))
    
    def forward(self, x):
        x = x.permute(0,2,1)
        x = F.relu(self.bn1(self.in_proj(x)))
        identity = x
        x = F.relu(self.bn2(self.conv1(x)))
        x = self.bn3(self.conv2(x))
        if identity.shape[1] != x.shape[1]:
            identity = F.conv1d(identity, torch.eye(256, 128, device=x.device).unsqueeze(2), padding=0)
        x = F.relu(x + identity)
        x = F.relu(self.bn4(self.conv3(x)))
        identity = x
        x = self.bn5(self.conv4(x))
        x = F.relu(x + identity)
        x = self.pool(x).view(x.shape[0], -1)
        x = self.proj(x)
        return F.normalize(x, dim=-1)

def seq_to_onehot(seq):
    import random
    arr = np.zeros((len(seq),4), dtype=np.float32)
    for i,ch in enumerate(seq):
        ch = ch.upper()
        if ch == "A": arr[i,0]=1
        elif ch=="C": arr[i,1]=1
        elif ch=="G": arr[i,2]=1
        elif ch=="T": arr[i,3]=1
        elif ch=="N": arr[i, random.randrange(4)] = 0.25
        else: arr[i, random.randrange(4)] = 1.0
    return arr

def load_reads(fastq_path, max_reads=None):
    """Load FASTQ reads (gzipped or plain)"""
    reads = []
    read_ids = []
    
    opener = gzip.open if fastq_path.endswith('.gz') else open
    mode = 'rt' if fastq_path.endswith('.gz') else 'r'
    
    print(f"  Loading reads from {os.path.basename(fastq_path)}...")
    with opener(fastq_path, mode) as f:
        for i, record in enumerate(SeqIO.parse(f, "fastq")):
            if max_reads and i >= max_reads:
                break
            
            seq = str(record.seq).upper()
            if len(seq) >= SEED_LEN:
                reads.append(seq[:SEED_LEN])  # use first SEED_LEN bases
                read_ids.append(record.id)
            
            if (i + 1) % 10000 == 0:
                print(f"    Loaded {i+1:,} reads...")
    
    return reads, read_ids

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Align HG002 reads to GRCh38')
    parser.add_argument('--reads', required=True, help='FASTQ file (can be .gz)')
    parser.add_argument('--max-reads', type=int, default=None, help='Max reads to process')
    parser.add_argument('--output', default='alignments.tsv', help='Output TSV file')
    parser.add_argument('--top-k', type=int, default=TOP_K, help='Top-K candidates')
    args = parser.parse_args()
    
    print("=" * 70)
    print("HG002 → GRCh38 Alignment Pipeline")
    print("=" * 70)
    
    # Load model
    print(f"\n[1/6] Loading encoder...")
    if not os.path.exists(ENCODER_PATH):
        print(f"ERROR: {ENCODER_PATH} not found")
        sys.exit(1)
    
    model = ImprovedNAL().to(DEVICE)
    model.load_state_dict(torch.load(ENCODER_PATH, map_location=DEVICE))
    model.eval()
    print(f"  ✓ Loaded on {DEVICE}")
    
    # Load index
    print(f"\n[2/6] Loading FAISS index...")
    if not os.path.exists(FAISS_INDEX):
        print(f"ERROR: {FAISS_INDEX} not found")
        print("Please run build_faiss_grch38.py first")
        sys.exit(1)
    
    index = faiss.read_index(FAISS_INDEX)
    print(f"  ✓ Loaded {index.ntotal:,} vectors")
    print(f"  ✓ nprobe = {index.nprobe}")
    
    # Load metadata
    print(f"\n[3/6] Loading reference metadata...")
    ref_pos = np.load(REF_POSITIONS)
    ref_chr = np.load(REF_CHROMOSOMES)
    print(f"  ✓ Positions: {len(ref_pos):,}")
    print(f"  ✓ Chromosomes: {len(np.unique(ref_chr))}")
    
    # Load reads
    print(f"\n[4/6] Loading query reads...")
    if not os.path.exists(args.reads):
        print(f"ERROR: {args.reads} not found")
        sys.exit(1)
    
    reads, read_ids = load_reads(args.reads, args.max_reads)
    print(f"  ✓ Loaded {len(reads):,} reads")
    
    # Encode reads
    print(f"\n[5/6] Encoding query reads...")
    all_vecs = []
    encode_start = time.time()
    
    with torch.no_grad():
        for i in range(0, len(reads), BATCH_SIZE):
            batch = reads[i:i+BATCH_SIZE]
            arr = np.stack([seq_to_onehot(s) for s in batch])
            tensor = torch.from_numpy(arr).to(DEVICE)
            vecs = model(tensor).cpu().numpy()
            all_vecs.append(vecs)
            
            if (i + BATCH_SIZE) % 5000 == 0 or i + BATCH_SIZE >= len(reads):
                progress = min(i + BATCH_SIZE, len(reads))
                elapsed = time.time() - encode_start
                rate = progress / elapsed
                print(f"    {progress:,}/{len(reads):,} ({progress/len(reads)*100:.1f}%) @ {rate:.0f} reads/s")
    
    query_vecs = np.vstack(all_vecs).astype('float32')
    encode_time = time.time() - encode_start
    print(f"  ✓ Encoded {len(query_vecs):,} queries in {encode_time:.1f}s ({len(query_vecs)/encode_time:.0f} reads/s)")
    
    # Search
    print(f"\n[6/6] Searching index (top-{args.top_k})...")
    search_start = time.time()
    D, I = index.search(query_vecs, args.top_k)
    search_time = time.time() - search_start
    print(f"  ✓ Search completed in {search_time:.2f}s ({len(query_vecs)/search_time:.0f} queries/s)")
    
    # Save results
    print(f"\nSaving alignments to {args.output}...")
    with open(args.output, 'w') as f:
        f.write("read_id\tchrom\tposition\tscore\trank\n")
        for qi in range(len(read_ids)):
            for rank in range(args.top_k):
                idx = I[qi, rank]
                if idx < 0 or idx >= len(ref_pos):
                    continue
                
                chr_id = ref_chr[idx]
                pos = ref_pos[idx]
                score = D[qi, rank]
                
                f.write(f"{read_ids[qi]}\t{chr_id}\t{pos}\t{score:.4f}\t{rank+1}\n")
    
    print(f"  ✓ Saved {len(read_ids) * args.top_k:,} alignments")
    
    # Statistics
    print(f"\n{'='*70}")
    print("Summary")
    print("=" * 70)
    print(f"Reads processed: {len(read_ids):,}")
    print(f"Total time: {time.time() - encode_start:.1f}s")
    print(f"  - Encoding: {encode_time:.1f}s ({len(query_vecs)/encode_time:.0f} reads/s)")
    print(f"  - Searching: {search_time:.2f}s ({len(query_vecs)/search_time:.0f} queries/s)")
    print(f"\nTop-1 alignment statistics:")
    
    # Chromosome distribution
    chr_counts = defaultdict(int)
    for qi in range(len(read_ids)):
        chr_id = ref_chr[I[qi, 0]]
        chr_counts[chr_id] += 1
    
    print(f"  Unique chromosomes hit: {len(chr_counts)}")
    top_chrs = sorted(chr_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"  Top chromosomes:")
    for chr_id, count in top_chrs:
        print(f"    Chr {chr_id}: {count:,} reads ({count/len(read_ids)*100:.1f}%)")
    
    print("=" * 70)

if __name__ == "__main__":
    main()
