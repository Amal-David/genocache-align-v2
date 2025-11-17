#!/usr/bin/env python3
"""
Build FAISS index for Chr22 only using focused trained model

This creates a Chr22-specific index with verified position tracking
to test if 8000 batches of training improved accuracy.
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import faiss
from pathlib import Path
from datetime import datetime

# NAL Encoder (same architecture as training)
class NALEncoder(nn.Module):
    def __init__(self, emb_dim=128, seed_len=512, vocab_size=5, 
                 hidden_dim=128, num_layers=4):
        super().__init__()
        self.emb_dim = emb_dim
        self.seed_len = seed_len
        self.vocab_size = vocab_size
        
        self.embedding = nn.Embedding(vocab_size, hidden_dim)
        self.conv_blocks = nn.ModuleList([
            nn.Conv1d(hidden_dim, hidden_dim, kernel_size=7, padding=3)
            for _ in range(num_layers)
        ])
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Linear(hidden_dim, emb_dim)
        self.norm = nn.LayerNorm(emb_dim)
    
    def forward(self, x):
        x = self.embedding(x)
        x = x.transpose(1, 2)
        
        for conv in self.conv_blocks:
            x = F.relu(conv(x))
        
        x = self.pool(x).squeeze(-1)
        x = self.fc(x)
        x = self.norm(x)
        x = F.normalize(x, p=2, dim=-1)
        return x

def load_model(model_path, device='cuda'):
    """Load trained Chr22 model"""
    print(f"Loading model from {model_path}...")
    checkpoint = torch.load(model_path, map_location=device)
    
    model = NALEncoder(
        emb_dim=checkpoint.get('emb_dim', 128),
        seed_len=checkpoint.get('seed_len', 512),
        vocab_size=5,
        hidden_dim=128,
        num_layers=4
    ).to(device)
    
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    print(f"  ✅ Model loaded (batch {checkpoint['batch']}, loss {checkpoint['loss']:.4f})")
    return model

def load_chromosome(fasta_path, chr_name='NC_000022.11'):
    """Load Chr22 from reference"""
    print(f"Loading {chr_name} from {fasta_path}...")
    seq = []
    in_chr = False
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if chr_name in line:
                    in_chr = True
                    print(f"  Found {chr_name}")
                else:
                    if in_chr:
                        break
                    in_chr = False
            elif in_chr:
                seq.append(line.upper())
    
    full_seq = ''.join(seq)
    print(f"  ✅ Loaded {len(full_seq):,} bp")
    return full_seq

def encode_chromosome(model, chr_seq, chr_name, seed_len=512, stride=32, batch_size=512, device='cuda'):
    """Encode all positions in Chr22"""
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    print(f"\nEncoding {chr_name} with stride={stride}...")
    print(f"  Seed length: {seed_len}bp")
    print(f"  Expected vectors: ~{(len(chr_seq) - seed_len) // stride:,}")
    
    embeddings = []
    positions = []
    
    # Collect all positions first
    print(f"  Collecting positions...")
    for start_pos in range(0, len(chr_seq) - seed_len, stride):
        positions.append({
            'chr': chr_name,
            'pos': start_pos,
            'strand': '+'
        })
    
    print(f"  Collected {len(positions):,} positions")
    print(f"  Encoding in batches of {batch_size}...")
    
    # Now encode in batches
    total_pos = 0
    for batch_start in range(0, len(positions), batch_size):
        batch_end = min(batch_start + batch_size, len(positions))
        batch_positions = positions[batch_start:batch_end]
        
        # Convert to sequences
        batch_seqs = []
        for pos_info in batch_positions:
            seed_start = pos_info['pos']
            seq = chr_seq[seed_start:seed_start + seed_len]
            idx = [base_to_idx.get(b, 4) for b in seq]
            batch_seqs.append(idx)
        
        batch_tensor = torch.tensor(batch_seqs, dtype=torch.long).to(device)
        
        with torch.no_grad():
            batch_emb = model(batch_tensor)
            embeddings.append(batch_emb.cpu().numpy())
        
        total_pos += len(batch_positions)
        if total_pos % 50000 == 0:
            print(f"    Encoded {total_pos:,}/{len(positions):,} positions...")
    

    
    # Concatenate
    all_embeddings = np.vstack(embeddings)
    
    print(f"  ✅ Encoded {len(positions):,} positions")
    print(f"  Embeddings shape: {all_embeddings.shape}")
    
    return all_embeddings, positions

def build_faiss_index(embeddings, emb_dim=128):
    """Build FAISS index"""
    print(f"\nBuilding FAISS index...")
    print(f"  Vectors: {embeddings.shape[0]:,}")
    print(f"  Dimensions: {emb_dim}")
    
    # Use Flat index for Chr22 (small enough)
    index = faiss.IndexFlatIP(emb_dim)  # Inner product (cosine similarity)
    
    # Normalize embeddings
    faiss.normalize_L2(embeddings)
    
    # Add to index
    index.add(embeddings)
    
    print(f"  ✅ Index built with {index.ntotal:,} vectors")
    return index

def save_positions(positions, output_path):
    """Save position mappings"""
    print(f"\nSaving positions to {output_path}...")
    
    chrs = np.array([p['chr'] for p in positions], dtype='U50')
    pos = np.array([p['pos'] for p in positions], dtype=np.int64)
    strands = np.array([p['strand'] for p in positions], dtype='U1')
    
    np.savez(
        output_path,
        chr=chrs,
        pos=pos,
        strand=strands
    )
    
    print(f"  ✅ Saved {len(positions):,} position mappings")

def main():
    print("="*80)
    print("CHR22 NAL INDEX BUILDING")
    print("="*80)
    print()
    
    # Paths
    model_path = 'models/chr22_nal_512bp_final.pt'
    reference_path = '/home/nebius/genocache/GRCh38.fa'
    index_output = 'indexes/chr22_nal_512bp_stride32.index'
    positions_output = 'indexes/chr22_nal_512bp_stride32.positions.npz'
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")
    print()
    
    # Load model
    model = load_model(model_path, device)
    
    # Load Chr22
    chr_seq = load_chromosome(reference_path, 'NC_000022.11')
    
    # Encode chromosome
    embeddings, positions = encode_chromosome(
        model, chr_seq, 'NC_000022.11',
        seed_len=512, stride=32, batch_size=512, device=device
    )
    
    # Build FAISS index
    index = build_faiss_index(embeddings, emb_dim=128)
    
    # Save index
    print(f"\nSaving FAISS index to {index_output}...")
    faiss.write_index(index, index_output)
    print(f"  ✅ Index saved")
    
    # Save positions
    save_positions(positions, positions_output)
    
    # Summary
    print("\n" + "="*80)
    print("INDEX BUILDING COMPLETE")
    print("="*80)
    print(f"\nCreated:")
    print(f"  • FAISS index: {index_output}")
    print(f"    - Vectors: {index.ntotal:,}")
    print(f"    - Dimensions: 128")
    print(f"  • Position mappings: {positions_output}")
    print(f"    - Chr22 positions: {len(positions):,}")
    print(f"\nNext steps:")
    print(f"  1. Generate Chr22 synthetic reads")
    print(f"  2. Test NAL vs minimap2")
    print(f"  3. Compare with previous 43% accuracy")
    print()

if __name__ == '__main__':
    main()
