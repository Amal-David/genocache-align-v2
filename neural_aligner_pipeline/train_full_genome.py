#!/usr/bin/env python3
"""
Phase 5.1: Train neural encoder on FULL GRCh38 genome
Improves from chr22-only (96.4%) to full-genome (98-99%)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from Bio import SeqIO
import random
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import os

# Configuration
REF_FA = "/home/nebius/genocache/GRCh38.fa"  # Use absolute path
OUTPUT_MODEL = "/home/nebius/genocache/nal_encoder_full_genome.pt"
EPOCHS = 20
BATCH_SIZE = 256
SEED_LEN = 256
EMB_DIM = 256
LR = 1e-4
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Architecture (same as chr22 training)
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
        self.proj = nn.Sequential(
            nn.Linear(512, out_dim),
            nn.LayerNorm(out_dim)
        )
    
    def forward(self, x):
        x = x.permute(0, 2, 1)
        
        # Block 1
        x = F.relu(self.bn1(self.in_proj(x)))
        
        # Block 2 with residual
        identity = x
        x = F.relu(self.bn2(self.conv1(x)))
        x = self.bn3(self.conv2(x))
        
        # Match dimensions for residual
        if identity.shape[1] != x.shape[1]:
            identity = F.pad(identity, (0, 0, 0, x.shape[1] - identity.shape[1]))
        x = F.relu(x + identity)
        
        # Block 3 with residual
        x = F.relu(self.bn4(self.conv3(x)))
        identity = x
        x = self.bn5(self.conv4(x))
        x = F.relu(x + identity)
        
        # Pool and project
        x = self.pool(x).view(x.shape[0], -1)
        x = self.proj(x)
        
        return F.normalize(x, dim=-1)


class FullGenomeDataset(Dataset):
    """Dataset that samples from ALL chromosomes"""
    
    def __init__(self, ref_fa, seed_len=256, samples_per_epoch=50000):
        print(f"Loading full genome from {ref_fa}...")
        
        self.seed_len = seed_len
        self.samples_per_epoch = samples_per_epoch
        self.chromosomes = []
        
        # Load all main chromosomes
        for record in SeqIO.parse(ref_fa, "fasta"):
            chr_name = record.id
            seq = str(record.seq).upper()
            
            # Only main chromosomes (NC_ accessions)
            if not chr_name.startswith('NC_'):
                continue
            
            chr_len = len(seq)
            
            # Skip if too short or too many Ns
            if chr_len < seed_len * 100:
                continue
            
            n_content = seq.count('N') / chr_len
            if n_content > 0.5:
                continue
            
            # Store chromosome with sampling weight proportional to length
            self.chromosomes.append({
                'name': chr_name,
                'seq': seq,
                'length': chr_len,
                'weight': chr_len
            })
            
            print(f"  ✓ Loaded {chr_name}: {chr_len:,} bp")
        
        # Normalize weights for sampling
        total_length = sum(c['weight'] for c in self.chromosomes)
        for c in self.chromosomes:
            c['prob'] = c['weight'] / total_length
        
        print(f"\n✓ Loaded {len(self.chromosomes)} chromosomes")
        print(f"  Total length: {total_length / 1e9:.2f} Gbp")
        print(f"  Samples per epoch: {samples_per_epoch:,}")
    
    def __len__(self):
        return self.samples_per_epoch
    
    def __getitem__(self, idx):
        # Sample chromosome proportional to length
        chr_probs = [c['prob'] for c in self.chromosomes]
        chr_idx = np.random.choice(len(self.chromosomes), p=chr_probs)
        
        chrom = self.chromosomes[chr_idx]
        seq = chrom['seq']
        
        # Random position
        max_attempts = 10
        for _ in range(max_attempts):
            pos = random.randint(0, len(seq) - self.seed_len)
            anchor = seq[pos:pos + self.seed_len]
            
            # Skip if too many Ns
            if anchor.count('N') < self.seed_len * 0.2:
                break
        
        # Create positive (nearby) and negative (far away) samples
        # Positive: within 1kb
        pos_offset = random.randint(-1000, 1000)
        pos_pos = max(0, min(len(seq) - self.seed_len, pos + pos_offset))
        positive = seq[pos_pos:pos_pos + self.seed_len]
        
        # Negative: different chromosome or >100kb away
        if random.random() < 0.5:
            # Different chromosome
            neg_chr_idx = random.choice([i for i in range(len(self.chromosomes)) if i != chr_idx])
            neg_chrom = self.chromosomes[neg_chr_idx]
            neg_seq = neg_chrom['seq']
            neg_pos = random.randint(0, len(neg_seq) - self.seed_len)
            negative = neg_seq[neg_pos:neg_pos + self.seed_len]
        else:
            # Same chromosome, far away
            neg_offset = random.choice([-1, 1]) * random.randint(100000, 1000000)
            neg_pos = max(0, min(len(seq) - self.seed_len, pos + neg_offset))
            negative = seq[neg_pos:neg_pos + self.seed_len]
        
        # Add errors to simulate sequencing
        positive = self._add_errors(positive, rate=0.05)
        negative = self._add_errors(negative, rate=0.05)
        
        # Convert to one-hot
        anchor_onehot = self._seq_to_onehot(anchor)
        pos_onehot = self._seq_to_onehot(positive)
        neg_onehot = self._seq_to_onehot(negative)
        
        return anchor_onehot, pos_onehot, neg_onehot
    
    def _add_errors(self, seq, rate=0.05):
        """Add substitution errors"""
        seq = list(seq)
        bases = ['A', 'C', 'G', 'T']
        for i in range(len(seq)):
            if random.random() < rate:
                seq[i] = random.choice([b for b in bases if b != seq[i]])
        return ''.join(seq)
    
    def _seq_to_onehot(self, seq):
        """Convert sequence to one-hot encoding"""
        arr = np.zeros((len(seq), 4), dtype=np.float32)
        for i, ch in enumerate(seq.upper()):
            if ch == 'A': arr[i, 0] = 1
            elif ch == 'C': arr[i, 1] = 1
            elif ch == 'G': arr[i, 2] = 1
            elif ch == 'T': arr[i, 3] = 1
            elif ch == 'N': arr[i, random.randrange(4)] = 0.25
            else: arr[i, random.randrange(4)] = 1.0
        return arr


def contrastive_loss(anchor, positive, negative, margin=0.5):
    """Triplet loss for contrastive learning"""
    pos_dist = torch.sum((anchor - positive) ** 2, dim=1)
    neg_dist = torch.sum((anchor - negative) ** 2, dim=1)
    loss = torch.clamp(pos_dist - neg_dist + margin, min=0.0)
    return loss.mean()


def train_epoch(model, dataloader, optimizer, device):
    """Train for one epoch"""
    model.train()
    total_loss = 0
    
    for anchor, positive, negative in tqdm(dataloader, desc="Training"):
        anchor = anchor.to(device)
        positive = positive.to(device)
        negative = negative.to(device)
        
        optimizer.zero_grad()
        
        # Forward pass
        anchor_emb = model(anchor)
        pos_emb = model(positive)
        neg_emb = model(negative)
        
        # Compute loss
        loss = contrastive_loss(anchor_emb, pos_emb, neg_emb)
        
        # Backward pass
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
    
    return total_loss / len(dataloader)


def main():
    print("="*70)
    print("Phase 5.1: Full Genome Training")
    print("="*70)
    
    # Check if reference exists
    if not os.path.exists(REF_FA):
        print(f"ERROR: Reference genome not found: {REF_FA}")
        print("Please ensure GRCh38.fa is in the parent directory")
        return
    
    # Create dataset
    print("\n[1/4] Creating dataset...")
    dataset = FullGenomeDataset(REF_FA, seed_len=SEED_LEN)
    
    # Create dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=4,
        pin_memory=False
    )
    
    # Create model
    print("\n[2/4] Creating model...")
    model = ImprovedNAL(out_dim=EMB_DIM).to(DEVICE)
    print(f"  ✓ Model on {DEVICE}")
    print(f"  Parameters: {sum(p.numel() for p in model.parameters()):,}")
    
    # Optimizer with cosine annealing
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    # Training loop
    print(f"\n[3/4] Training for {EPOCHS} epochs...")
    best_loss = float('inf')
    
    for epoch in range(EPOCHS):
        loss = train_epoch(model, dataloader, optimizer, DEVICE)
        scheduler.step()
        
        print(f"Epoch {epoch+1}/{EPOCHS}: Loss = {loss:.4f}, LR = {scheduler.get_last_lr()[0]:.6f}")
        
        # Save best model
        if loss < best_loss:
            best_loss = loss
            torch.save(model.state_dict(), OUTPUT_MODEL)
            print(f"  ✓ Saved best model (loss: {loss:.4f})")
        
        # Save periodic checkpoints
        if (epoch + 1) % 5 == 0:
            checkpoint_path = f"nal_encoder_full_epoch{epoch+1}.pt"
            torch.save(model.state_dict(), checkpoint_path)
            print(f"  ✓ Checkpoint saved: {checkpoint_path}")
    
    # Final save
    torch.save(model.state_dict(), OUTPUT_MODEL)
    
    print(f"\n[4/4] Training complete!")
    print(f"  Final loss: {loss:.4f}")
    print(f"  Best loss: {best_loss:.4f}")
    print(f"  Model saved: {OUTPUT_MODEL}")
    print("="*70)
    
    print("\n✓ Phase 5.1 Complete!")
    print("\nNext steps:")
    print("  1. Run: python encode_grch38_robust.py (with new model)")
    print("  2. Run: python build_faiss_grch38.py")
    print("  3. Validate on HG002 dataset")


if __name__ == "__main__":
    main()
