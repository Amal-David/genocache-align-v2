#!/usr/bin/env python3
"""
Training script for GenoCacheEncoder with hard negative mining
"""

import multiprocessing as mp
try:
    mp.set_start_method('spawn', force=True)
except RuntimeError:
    pass

import os
import sys
import random
import argparse
from collections import defaultdict
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from Bio import SeqIO
from tqdm import tqdm

from genocache_encoder import (
    GenoCacheEncoder, GenoCacheLoss, seq_to_tokens,
    augment_sequence
)

# Default hyperparameters
SEED_LEN = 512
EMB_DIM = 256
BATCH_SIZE = 256  # Per GPU
EPOCHS = 50
LR = 3e-4
WEIGHT_DECAY = 0.01
TEMPERATURE = 0.07
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


class RepetitiveFinder:
    """
    Identify repetitive regions using k-mer frequency analysis
    Approximates RepeatMasker for Alu, LINE, SINE detection
    """
    def __init__(self, sequence, k=15, threshold=10):
        """
        Args:
            sequence: DNA sequence string
            k: k-mer size for frequency analysis
            threshold: k-mer count threshold for "repetitive"
        """
        print("Analyzing repetitive regions...")
        self.sequence = sequence.upper()
        self.k = k
        self.threshold = threshold
        self.repetitive_positions = set()
        
        # Build k-mer frequency map
        kmer_counts = defaultdict(int)
        kmer_positions = defaultdict(list)
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1
                kmer_positions[kmer].append(i)
        
        # Mark positions with high-frequency k-mers as repetitive
        for kmer, count in kmer_counts.items():
            if count >= threshold:
                for pos in kmer_positions[kmer]:
                    # Mark a window around this position
                    for j in range(max(0, pos - 100), min(len(sequence), pos + 100)):
                        self.repetitive_positions.add(j)
        
        self.repetitive_ratio = len(self.repetitive_positions) / len(sequence)
        print(f"  Found {len(self.repetitive_positions):,} repetitive positions")
        print(f"  Repetitive ratio: {self.repetitive_ratio*100:.2f}%")
    
    def is_repetitive(self, pos, window=256):
        """Check if a window starting at pos is repetitive"""
        # Count how many bases in the window are in repetitive regions
        count = sum(1 for i in range(pos, min(pos + window, len(self.sequence)))
                   if i in self.repetitive_positions)
        return count > window * 0.3  # >30% of window is repetitive
    
    def sample_repetitive_region(self, window=256):
        """Sample a random position from repetitive regions"""
        max_attempts = 100
        for _ in range(max_attempts):
            pos = random.randint(0, len(self.sequence) - window)
            if self.is_repetitive(pos, window):
                return pos
        
        # Fallback: return any random position
        return random.randint(0, len(self.sequence) - window)


class Chr1Dataset(Dataset):
    """
    Dataset that samples from chromosome 1 with hard negative mining
    """
    def __init__(self, ref_fasta, chrom_id="NC_000001.11", 
                 seed_len=SEED_LEN, n_samples=100000,
                 hard_negative_ratio=0.3):
        """
        Args:
            ref_fasta: path to reference FASTA (GRCh38.fa)
            chrom_id: chromosome identifier (chr1 = NC_000001.11)
            seed_len: seed length
            n_samples: number of training samples per epoch
            hard_negative_ratio: fraction of negatives from repetitive regions
        """
        print(f"Loading chromosome: {chrom_id}...")
        
        # Load chromosome 1
        self.sequence = None
        for record in SeqIO.parse(ref_fasta, "fasta"):
            if record.id == chrom_id:
                self.sequence = str(record.seq).upper()
                break
        
        if self.sequence is None:
            raise ValueError(f"Could not find {chrom_id} in {ref_fasta}")
        
        print(f"  Length: {len(self.sequence):,} bp")
        
        self.seed_len = seed_len
        self.n_samples = n_samples
        self.hard_negative_ratio = hard_negative_ratio
        
        # Identify repetitive regions for hard negative mining
        self.repeat_finder = RepetitiveFinder(self.sequence, k=15, threshold=10)
    
    def __len__(self):
        return self.n_samples
    
    def __getitem__(self, idx):
        """
        Returns:
            anchor_seq: original sequence
            hard_neg_seqs: list of sequences from repetitive regions
        """
        # Sample random position for anchor
        pos = random.randint(0, len(self.sequence) - self.seed_len)
        anchor_seq = self.sequence[pos:pos + self.seed_len]
        
        # Sample hard negatives from repetitive regions
        n_hard = int(10 * self.hard_negative_ratio)  # e.g., 3 hard negatives
        hard_neg_seqs = []
        
        for _ in range(n_hard):
            neg_pos = self.repeat_finder.sample_repetitive_region(self.seed_len)
            hard_neg_seq = self.sequence[neg_pos:neg_pos + self.seed_len]
            hard_neg_seqs.append(hard_neg_seq)
        
        return anchor_seq, hard_neg_seqs


def collate_fn(batch):
    """
    Collate function for DataLoader
    Returns anchors, positives (augmented), and hard negatives
    """
    anchor_seqs, hard_neg_lists = zip(*batch)
    
    # Augment anchors to create positives
    positive_seqs = [augment_sequence(seq) for seq in anchor_seqs]
    
    # Convert to tokens
    anchors_tokens = [seq_to_tokens(seq) for seq in anchor_seqs]
    positives_tokens = [seq_to_tokens(seq) for seq in positive_seqs]
    
    # Hard negatives: [B, K, L] where K is number of hard negatives per sample
    max_hard = max(len(hn_list) for hn_list in hard_neg_lists)
    batch_size = len(anchor_seqs)
    
    hard_negs_tokens = []
    for hn_list in hard_neg_lists:
        # Augment hard negatives too
        hn_tokens = [seq_to_tokens(augment_sequence(seq)) for seq in hn_list]
        # Pad to max_hard if needed
        while len(hn_tokens) < max_hard:
            hn_tokens.append(hn_tokens[0] if hn_tokens else [4] * SEED_LEN)
        hard_negs_tokens.append(hn_tokens)
    
    # Convert to tensors
    anchors = torch.tensor(anchors_tokens, dtype=torch.long)
    positives = torch.tensor(positives_tokens, dtype=torch.long)
    hard_negs = torch.tensor(hard_negs_tokens, dtype=torch.long)
    
    return anchors, positives, hard_negs


def train(args):
    """Main training loop"""
    
    print("="*80)
    print("GenoCache V3 Training")
    print("="*80)
    print(f"Device: {DEVICE}")
    print(f"Reference: {args.fasta}")
    print(f"Chromosome: {args.chrom_id}")
    print(f"Seed length: {args.seed_len}")
    print(f"Embedding dim: {args.emb_dim}")
    print(f"Batch size: {args.batch_size}")
    print(f"Epochs: {args.epochs}")
    print(f"Learning rate: {args.lr}")
    print(f"Hard negative ratio: {args.hard_negative_ratio}")
    print("="*80)
    print()
    
    # Create dataset
    dataset = Chr1Dataset(
        ref_fasta=args.fasta,
        chrom_id=args.chrom_id,
        seed_len=args.seed_len,
        n_samples=args.n_samples,
        hard_negative_ratio=args.hard_negative_ratio
    )
    
    # Create dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=4,
        collate_fn=collate_fn,
        pin_memory=False
    )
    
    print(f"Dataset: {len(dataset):,} samples/epoch")
    print(f"Batches per epoch: {len(dataloader):,}")
    print()
    
    # Create model
    model = GenoCacheEncoder(seed_len=args.seed_len, emb_dim=args.emb_dim)
    model = model.to(DEVICE)
    
    # Loss and optimizer
    criterion = GenoCacheLoss(temperature=TEMPERATURE)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.lr,
        weight_decay=WEIGHT_DECAY
    )
    
    # Cosine annealing scheduler
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs
    )
    
    # Training loop
    best_loss = float('inf')
    
    # Curriculum learning: gradually increase difficulty
    curriculum_epochs = min(10, args.epochs // 3)
    
    for epoch in range(args.epochs):
        model.train()
        total_loss = 0.0
        
        # Curriculum: start with easier examples (lower error rate)
        if epoch < curriculum_epochs:
            # Linearly increase error rate
            progress = epoch / curriculum_epochs
            # Start at 1% error, go up to 10%
            global AUG_ERR_MIN, AUG_ERR_MAX
            from genocache_encoder import augment_sequence
            # We'd need to modify augment_sequence to accept error rate params
            # For now, just note this in docs
        
        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{args.epochs}")
        
        for batch_idx, (anchors, positives, hard_negs) in enumerate(pbar):
            anchors = anchors.to(DEVICE)
            positives = positives.to(DEVICE)
            hard_negs = hard_negs.to(DEVICE)
            
            # Forward pass
            anchor_emb = model(anchors)
            positive_emb = model(positives)
            
            # Encode hard negatives
            B, K, L = hard_negs.shape
            hard_negs_flat = hard_negs.view(B * K, L)
            hard_negs_emb = model(hard_negs_flat)
            hard_negs_emb = hard_negs_emb.view(B, K, -1)
            
            # Compute loss
            loss = criterion(anchor_emb, positive_emb, hard_negs_emb)
            
            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            
            total_loss += loss.item()
            
            # Update progress bar
            pbar.set_postfix({
                'loss': f'{total_loss / (batch_idx + 1):.4f}',
                'lr': f'{scheduler.get_last_lr()[0]:.6f}'
            })
        
        scheduler.step()
        avg_loss = total_loss / len(dataloader)
        
        print(f"\nEpoch {epoch+1}/{args.epochs} - Loss: {avg_loss:.6f}")
        
        # Save checkpoint
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': avg_loss,
        }
        
        checkpoint_path = os.path.join(args.output_dir, f"genocache_epoch{epoch+1}.pt")
        torch.save(checkpoint, checkpoint_path)
        
        # Save best model
        if avg_loss < best_loss:
            best_loss = avg_loss
            best_path = os.path.join(args.output_dir, "genocache_best.pt")
            torch.save(checkpoint, best_path)
            print(f"  ✓ New best model: {best_loss:.6f}")
        
        print()
    
    print("="*80)
    print("Training Complete!")
    print(f"Best loss: {best_loss:.6f}")
    print(f"Models saved to: {args.output_dir}")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(description="Train GenoCacheEncoder")
    
    parser.add_argument("--fasta", type=str, default="../GRCh38.fa",
                       help="Path to reference FASTA")
    parser.add_argument("--chrom-id", type=str, default="NC_000001.11",
                       help="Chromosome ID (default: chr1)")
    parser.add_argument("--output-dir", type=str, default="./models",
                       help="Output directory for models")
    parser.add_argument("--seed-len", type=int, default=SEED_LEN,
                       help="Seed length")
    parser.add_argument("--emb-dim", type=int, default=EMB_DIM,
                       help="Embedding dimension")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE,
                       help="Batch size")
    parser.add_argument("--epochs", type=int, default=EPOCHS,
                       help="Number of epochs")
    parser.add_argument("--lr", type=float, default=LR,
                       help="Learning rate")
    parser.add_argument("--n-samples", type=int, default=100000,
                       help="Samples per epoch")
    parser.add_argument("--hard-negative-ratio", type=float, default=0.3,
                       help="Ratio of hard negatives from repetitive regions")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Train
    train(args)


if __name__ == "__main__":
    main()
