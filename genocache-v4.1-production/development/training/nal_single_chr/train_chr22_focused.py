#!/usr/bin/env python3
"""
Focused NAL Training - Single Chromosome with Curriculum Learning

Key improvements:
1. Train ONLY on Chr22 (controlled subset)
2. 4000-8000 batches (vs 62 previously)
3. Curriculum learning: 5% → 15% error rate
4. Adaptive seeding: 256bp + 512bp seeds
5. Careful validation at each step

Based on original NAL training but focused and improved.
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from datetime import datetime
import random

# Encoder model (same as before)
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
    print(f"  Loaded {len(full_seq):,} bp")
    return full_seq

def add_errors(seq, error_rate, error_types=['sub', 'ins', 'del']):
    """Add realistic sequencing errors"""
    result = []
    i = 0
    
    while i < len(seq):
        if random.random() < error_rate:
            error_type = random.choice(error_types)
            
            if error_type == 'sub':
                result.append(random.choice('ACGT'))
                i += 1
            elif error_type == 'ins':
                result.append(random.choice('ACGT'))
            else:  # del
                i += 1
        else:
            result.append(seq[i])
            i += 1
    
    return ''.join(result)

def generate_training_batch(chr_seq, batch_size, seed_len, error_rate, device):
    """Generate training batch with curriculum learning"""
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    anchors = []
    positives = []
    
    max_start = len(chr_seq) - seed_len - 1
    
    for _ in range(batch_size):
        # Random position
        start = random.randint(0, max_start)
        anchor_seq = chr_seq[start:start + seed_len]
        
        # Convert to indices
        anchor_idx = [base_to_idx.get(b, 4) for b in anchor_seq]
        
        # Add errors for positive
        pos_seq = add_errors(anchor_seq, error_rate)
        pos_idx = [base_to_idx.get(b, 4) for b in pos_seq[:seed_len]]
        
        # Pad if needed
        while len(pos_idx) < seed_len:
            pos_idx.append(4)
        pos_idx = pos_idx[:seed_len]
        
        anchors.append(anchor_idx)
        positives.append(pos_idx)
    
    return (
        torch.tensor(anchors, dtype=torch.long).to(device),
        torch.tensor(positives, dtype=torch.long).to(device)
    )

def info_nce_loss(anchor_emb, pos_emb, temperature=0.07):
    """InfoNCE contrastive loss"""
    # Normalize
    anchor_emb = F.normalize(anchor_emb, p=2, dim=-1)
    pos_emb = F.normalize(pos_emb, p=2, dim=-1)
    
    # Compute similarities
    similarity = torch.matmul(anchor_emb, pos_emb.T) / temperature
    
    # Labels (diagonal)
    labels = torch.arange(similarity.size(0)).to(similarity.device)
    
    # Cross entropy
    loss = F.cross_entropy(similarity, labels)
    
    # Compute accuracy
    pred = similarity.argmax(dim=1)
    acc = (pred == labels).float().mean()
    
    return loss, acc

def train_focused_chr22(
    reference_path,
    output_dir,
    num_batches=4000,
    batch_size=32,
    emb_dim=128,
    seed_lens=[256, 512],
    device='cuda'
):
    """
    Train NAL on Chr22 only with curriculum learning
    """
    print("="*80)
    print("FOCUSED CHR22 NAL TRAINING")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Chromosome: Chr22 (NC_000022.11)")
    print(f"  Batches: {num_batches}")
    print(f"  Batch size: {batch_size}")
    print(f"  Embedding dim: {emb_dim}")
    print(f"  Seed lengths: {seed_lens}")
    print(f"  Device: {device}")
    print()
    
    # Load Chr22
    chr_seq = load_chromosome(reference_path, 'NC_000022.11')
    
    if len(chr_seq) < 1_000_000:
        print("❌ Chr22 too short or not found!")
        return None
    
    # Create models for each seed length
    models = {}
    optimizers = {}
    
    for seed_len in seed_lens:
        print(f"\nInitializing model for {seed_len}bp seeds...")
        model = NALEncoder(
            emb_dim=emb_dim,
            seed_len=seed_len,
            vocab_size=5,
            hidden_dim=128,
            num_layers=4
        ).to(device)
        
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=1e-4,
            weight_decay=0.01
        )
        
        models[seed_len] = model
        optimizers[seed_len] = optimizer
    
    # Training with curriculum learning
    print("\n" + "="*80)
    print("TRAINING WITH CURRICULUM LEARNING")
    print("="*80)
    
    # Curriculum: gradually increase error rate (extended for 8000 batches)
    curriculum_schedule = [
        (0, 1000, 0.05),       # Batches 0-1000: 5% error
        (1000, 2500, 0.08),    # Batches 1000-2500: 8% error
        (2500, 5000, 0.12),    # Batches 2500-5000: 12% error
        (5000, num_batches, 0.15)  # Batches 5000+: 15% error
    ]
    
    def get_error_rate(batch_idx):
        for start, end, rate in curriculum_schedule:
            if start <= batch_idx < end:
                return rate
        return 0.15
    
    best_loss = float('inf')
    
    for batch_idx in range(num_batches):
        # Curriculum learning: increase error rate
        error_rate = get_error_rate(batch_idx)
        
        # Alternate between seed lengths
        seed_len = random.choice(seed_lens)
        model = models[seed_len]
        optimizer = optimizers[seed_len]
        
        # Generate batch
        anchor, positive = generate_training_batch(
            chr_seq, batch_size, seed_len, error_rate, device
        )
        
        # Forward pass
        model.train()
        anchor_emb = model(anchor)
        pos_emb = model(positive)
        
        # Loss
        loss, acc = info_nce_loss(anchor_emb, pos_emb)
        
        # Backward
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        
        # Logging
        if (batch_idx + 1) % 100 == 0:
            print(f"Batch {batch_idx+1}/{num_batches} | "
                  f"Seed: {seed_len}bp | Error: {error_rate:.1%} | "
                  f"Loss: {loss.item():.4f} | Acc: {acc.item():.1%}")
        
        # Save checkpoint
        if (batch_idx + 1) % 1000 == 0:
            if loss.item() < best_loss:
                best_loss = loss.item()
                
                for sl, model in models.items():
                    checkpoint_path = Path(output_dir) / f'models/chr22_nal_{sl}bp_batch{batch_idx+1}.pt'
                    torch.save({
                        'batch': batch_idx + 1,
                        'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': optimizers[sl].state_dict(),
                        'loss': loss.item(),
                        'acc': acc.item(),
                        'seed_len': sl,
                        'emb_dim': emb_dim
                    }, checkpoint_path)
                
                print(f"  ✅ Checkpoint saved at batch {batch_idx+1}")
    
    # Save final models
    print("\n" + "="*80)
    print("SAVING FINAL MODELS")
    print("="*80)
    
    for seed_len, model in models.items():
        final_path = Path(output_dir) / f'models/chr22_nal_{seed_len}bp_final.pt'
        torch.save({
            'batch': num_batches,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizers[seed_len].state_dict(),
            'loss': loss.item(),
            'acc': acc.item(),
            'seed_len': seed_len,
            'emb_dim': emb_dim,
            'chr': 'NC_000022.11'
        }, final_path)
        print(f"  ✅ Saved {seed_len}bp model: {final_path}")
    
    print("\n✅ Training complete!")
    return models

def main():
    reference_path = '/home/nebius/genocache/GRCh38.fa'
    output_dir = '/home/nebius/genocache/genocache-v4.1-production/development/training/nal_single_chr'
    
    # Train
    models = train_focused_chr22(
        reference_path=reference_path,
        output_dir=output_dir,
        num_batches=8000,  # Extended to 8000 for single chromosome
        batch_size=32,
        emb_dim=128,
        seed_lens=[256, 512],  # Adaptive seeding
        device='cuda'
    )
    
    if models:
        print("\n" + "="*80)
        print("NEXT STEPS:")
        print("="*80)
        print("1. Build index for Chr22 only")
        print("2. Test on Chr22 reads")
        print("3. Compare with minimap2")
        print("4. Check false positive rate")
        print()

if __name__ == '__main__':
    main()
