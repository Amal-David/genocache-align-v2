#!/usr/bin/env python3
"""
GenoCache-Align Training Pipeline (Single GPU Optimized)
Adapted from Claude's multi-GPU implementation for single H100
"""

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
from pathlib import Path
from tqdm import tqdm
import argparse
from datetime import datetime
import json
import os
import sys

from genocache_encoder import GenoCacheEncoder, GenoCacheLoss, augment_sequence, augment_rc


class GenomicSeedDataset(Dataset):
    """
    Enhanced dataset with hard negative mining from repetitive regions
    """
    
    def __init__(
        self,
        fasta_path,
        seed_len=512,
        samples_per_epoch=100000,
        hard_negative_ratio=0.3,
        curriculum_epoch=0,
        max_curriculum_epochs=10
    ):
        self.seed_len = seed_len
        self.samples_per_epoch = samples_per_epoch
        self.hard_negative_ratio = hard_negative_ratio
        self.curriculum_epoch = curriculum_epoch
        self.max_curriculum_epochs = max_curriculum_epochs
        
        # Load reference genome
        from Bio import SeqIO
        self.sequences = {}
        self.repetitive_regions = {}
        
        print(f"\n{'='*80}")
        print("Loading Reference Genome")
        print(f"{'='*80}")
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            # Filter out sequences with too many Ns
            if seq.count('N') / len(seq) < 0.1:
                self.sequences[record.id] = seq
                print(f"  Loading {record.id}... ", end='', flush=True)
                self.repetitive_regions[record.id] = self._find_repetitive_regions(seq)
                print(f"✓ ({len(seq):,} bp, {len(self.repetitive_regions[record.id])} repetitive regions)")
        
        self.chrom_names = list(self.sequences.keys())
        self.chrom_lengths = [len(self.sequences[chrom]) for chrom in self.chrom_names]
        self.total_length = sum(self.chrom_lengths)
        
        print(f"\n✅ Loaded {len(self.sequences)} chromosomes")
        print(f"   Total length: {self.total_length:,} bp")
        print(f"   Repetitive regions: {sum(len(r) for r in self.repetitive_regions.values())} identified")
        print(f"   Hard negative ratio: {hard_negative_ratio:.2%}")
        print(f"   Curriculum epoch: {curriculum_epoch}/{max_curriculum_epochs}")
        print(f"{'='*80}\n")
    
    def _find_repetitive_regions(self, seq, kmer_size=15, threshold=10):
        """Identify repetitive regions by counting k-mer occurrences"""
        from collections import defaultdict
        
        kmer_counts = defaultdict(int)
        kmer_positions = defaultdict(list)
        
        # Count k-mers (sample every 10bp to speed up)
        for i in range(0, len(seq) - kmer_size + 1, 10):
            kmer = seq[i:i+kmer_size]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1
                kmer_positions[kmer].append(i)
        
        # Find repetitive k-mers
        repetitive_regions = []
        for kmer, count in kmer_counts.items():
            if count >= threshold:
                for pos in kmer_positions[kmer]:
                    start = max(0, pos - self.seed_len // 2)
                    end = min(len(seq), pos + kmer_size + self.seed_len // 2)
                    repetitive_regions.append((start, end))
        
        # Merge overlapping regions
        if repetitive_regions:
            repetitive_regions.sort()
            merged = [repetitive_regions[0]]
            for start, end in repetitive_regions[1:]:
                if start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            return merged
        return []
    
    def __len__(self):
        return self.samples_per_epoch
    
    def _get_curriculum_error_rate(self):
        """Curriculum learning: start with low error rates, gradually increase"""
        if self.curriculum_epoch >= self.max_curriculum_epochs:
            return np.random.uniform(0.01, 0.10)
        else:
            progress = self.curriculum_epoch / self.max_curriculum_epochs
            min_error = 0.01
            max_error = 0.01 + (0.10 - 0.01) * progress
            return np.random.uniform(min_error, max_error)
    
    def _sample_hard_negative(self, chrom, exclude_start, exclude_end):
        """Sample a hard negative from repetitive regions"""
        repetitive = self.repetitive_regions[chrom]
        
        if not repetitive:
            return self._sample_random_seed(chrom)
        
        # Try to sample from repetitive region, fallback to random if fails
        max_attempts = 10
        for attempt in range(max_attempts):
            region = repetitive[np.random.randint(0, len(repetitive))]
            region_start, region_end = region
            
            # Ensure region is large enough
            if region_end - region_start < self.seed_len:
                continue
            
            # Try to find valid position avoiding anchor
            valid_ranges = []
            
            # Before anchor (if far enough)
            if region_start < exclude_start - self.seed_len * 2:
                end_pos = min(exclude_start - self.seed_len * 2, region_end - self.seed_len)
                if end_pos > region_start:
                    valid_ranges.append((region_start, end_pos))
            
            # After anchor (if far enough)
            if region_end > exclude_end + self.seed_len * 2:
                start_pos = max(exclude_end + self.seed_len * 2, region_start)
                end_pos = region_end - self.seed_len
                if end_pos > start_pos:
                    valid_ranges.append((start_pos, end_pos))
            
            # If no valid ranges avoiding anchor, use whole region
            if not valid_ranges:
                if region_end - self.seed_len > region_start:
                    valid_ranges.append((region_start, region_end - self.seed_len))
            
            # Sample from valid ranges
            if valid_ranges:
                range_idx = np.random.randint(0, len(valid_ranges))
                start, end = valid_ranges[range_idx]
                if end > start:
                    pos = np.random.randint(start, end)
                    seq = self.sequences[chrom][pos:pos + self.seed_len]
                    
                    # Check quality
                    if seq.count('N') <= self.seed_len * 0.1:
                        return seq
        
        # Fallback to random seed
        return self._sample_random_seed(chrom)
    
    def _sample_random_seed(self, chrom):
        """Sample a random seed from chromosome"""
        chrom_seq = self.sequences[chrom]
        max_pos = len(chrom_seq) - self.seed_len
        if max_pos <= 0:
            return None
        pos = np.random.randint(0, max_pos)
        return chrom_seq[pos:pos + self.seed_len]
    
    def __getitem__(self, idx):
        # Sample random chromosome weighted by length
        chrom_idx = np.random.choice(
            len(self.chrom_names),
            p=np.array(self.chrom_lengths) / self.total_length
        )
        chrom = self.chrom_names[chrom_idx]
        chrom_seq = self.sequences[chrom]
        
        # Sample random position
        max_pos = len(chrom_seq) - self.seed_len
        if max_pos <= 0:
            return self.__getitem__(idx)
        
        pos = np.random.randint(0, max_pos)
        seed = chrom_seq[pos:pos + self.seed_len]
        
        # Skip if too many Ns
        if seed.count('N') > self.seed_len * 0.1:
            return self.__getitem__(idx)
        
        # Generate augmented version (positive pair)
        error_rate = self._get_curriculum_error_rate()
        augmented = augment_sequence(seed, error_rate=error_rate, seed_len=self.seed_len)
        
        # RC augmentation
        seed = augment_rc(seed, prob=0.5)
        augmented = augment_rc(augmented, prob=0.5)
        
        # Sample hard negatives
        hard_negatives = []
        if np.random.random() < self.hard_negative_ratio:
            num_hard_negs = np.random.randint(2, 5)
            for _ in range(num_hard_negs):
                hard_neg = self._sample_hard_negative(chrom, pos, pos + self.seed_len)
                if hard_neg is not None:
                    hard_neg = augment_rc(hard_neg, prob=0.5)
                    hard_negatives.append(hard_neg)
        
        return seed, augmented, hard_negatives


def collate_fn(batch, seed_len=512):
    """Convert DNA strings to tensors with variable hard negatives"""
    seeds, augmented, hard_negatives_list = zip(*batch)
    
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def str_to_tensor(seq):
        """Convert sequence to tensor, ensuring fixed length"""
        # Pad or trim to seed_len
        if len(seq) < seed_len:
            seq = seq + 'N' * (seed_len - len(seq))
        elif len(seq) > seed_len:
            seq = seq[:seed_len]
        return torch.tensor([mapping.get(base, 4) for base in seq.upper()], dtype=torch.long)
    
    seed_tensors = torch.stack([str_to_tensor(s) for s in seeds])
    aug_tensors = torch.stack([str_to_tensor(a) for a in augmented])
    
    # Handle hard negatives
    max_hard_negs = max(len(hn) for hn in hard_negatives_list)
    if max_hard_negs > 0:
        hard_neg_tensors = []
        for hard_negs in hard_negatives_list:
            if hard_negs:
                tensors = [str_to_tensor(hn) for hn in hard_negs]
                while len(tensors) < max_hard_negs:
                    tensors.append(tensors[0].clone())
                hard_neg_tensors.append(torch.stack(tensors))
            else:
                # Create dummy hard negatives
                dummy = seed_tensors[:1].repeat(max_hard_negs, 1)
                hard_neg_tensors.append(dummy)
        hard_neg_tensors = torch.stack(hard_neg_tensors)
    else:
        hard_neg_tensors = None
    
    return seed_tensors, aug_tensors, hard_neg_tensors


def train_epoch(model, dataloader, criterion, optimizer, device, epoch):
    """Train for one epoch"""
    model.train()
    total_loss = 0
    num_batches = 0
    
    pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}")
    
    for batch_idx, (seeds, augmented, hard_negatives) in enumerate(pbar):
        seeds = seeds.to(device)
        augmented = augmented.to(device)
        if hard_negatives is not None:
            hard_negatives = hard_negatives.to(device)
        
        optimizer.zero_grad()
        
        # Forward pass
        embeddings = model(seeds, use_contrastive_head=True)
        aug_embeddings = model(augmented, use_contrastive_head=True)
        
        # Encode hard negatives
        hard_neg_embeddings = None
        if hard_negatives is not None:
            B, K, L = hard_negatives.shape
            hard_negatives_flat = hard_negatives.view(B * K, L)
            hard_neg_embeddings_flat = model(hard_negatives_flat, use_contrastive_head=True)
            hard_neg_embeddings = hard_neg_embeddings_flat.view(B, K, -1)
        
        # Compute loss
        loss = criterion(embeddings, aug_embeddings, hard_neg_embeddings)
        
        # Backward pass
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        total_loss += loss.item()
        num_batches += 1
        
        pbar.set_postfix({
            'loss': f'{loss.item():.4f}',
            'avg_loss': f'{total_loss/num_batches:.4f}'
        })
    
    avg_loss = total_loss / num_batches
    return avg_loss


def main():
    parser = argparse.ArgumentParser(description='Train GenoCacheEncoder (Single GPU)')
    
    # Data
    parser.add_argument('--fasta', type=str, required=True,
                       help='Path to reference genome FASTA')
    parser.add_argument('--output-dir', type=str, default='./models',
                       help='Output directory')
    
    # Model
    parser.add_argument('--seed-len', type=int, default=512,
                       help='Seed length')
    parser.add_argument('--emb-dim', type=int, default=256,
                       help='Embedding dimension')
    
    # Training
    parser.add_argument('--epochs', type=int, default=50,
                       help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=256,
                       help='Batch size')
    parser.add_argument('--samples-per-epoch', type=int, default=100000,
                       help='Samples per epoch')
    parser.add_argument('--learning-rate', type=float, default=1e-3,
                       help='Initial learning rate')
    parser.add_argument('--weight-decay', type=float, default=0.1,
                       help='Weight decay')
    parser.add_argument('--temperature', type=float, default=0.07,
                       help='Temperature for InfoNCE loss')
    parser.add_argument('--patience', type=int, default=12,
                       help='Early stopping patience')
    
    # Advanced
    parser.add_argument('--hard-negative-ratio', type=float, default=0.3,
                       help='Ratio of samples with hard negatives')
    parser.add_argument('--hard-negative-weight', type=float, default=2.0,
                       help='Weight for hard negative loss')
    parser.add_argument('--max-curriculum-epochs', type=int, default=10,
                       help='Number of curriculum learning epochs')
    
    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n{'='*80}")
    print("GENOCACHE-ALIGN ENCODER TRAINING (Single GPU)")
    print(f"{'='*80}")
    print(f"Device: {device}")
    print(f"Seed length: {args.seed_len}")
    print(f"Embedding dim: {args.emb_dim}")
    print(f"Batch size: {args.batch_size}")
    print(f"Epochs: {args.epochs}")
    print(f"Learning rate: {args.learning_rate}")
    print(f"Hard negative mining: {args.hard_negative_ratio:.2%}")
    print(f"Curriculum learning: {args.max_curriculum_epochs} epochs")
    print(f"{'='*80}\n")
    
    # Create output directories
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / 'checkpoints'
    checkpoint_dir.mkdir(exist_ok=True)
    
    # Create model
    model = GenoCacheEncoder(
        emb_dim=args.emb_dim,
        seed_len=args.seed_len
    ).to(device)
    
    model_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model parameters: {model_params:,}\n")
    
    # Loss and optimizer
    criterion = GenoCacheLoss(
        temperature=args.temperature,
        hard_negative_weight=args.hard_negative_weight
    )
    
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay
    )
    
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=3
    )
    
    # Training loop
    best_loss = float('inf')
    no_improve_epochs = 0
    
    training_log = {
        'model': 'GenoCacheEncoder',
        'config': vars(args),
        'start_time': datetime.now().isoformat(),
        'epochs': []
    }
    
    for epoch in range(args.epochs):
        curriculum_epoch = min(epoch, args.max_curriculum_epochs)
        
        # Create dataset
        train_dataset = GenomicSeedDataset(
            args.fasta,
            seed_len=args.seed_len,
            samples_per_epoch=args.samples_per_epoch,
            hard_negative_ratio=args.hard_negative_ratio,
            curriculum_epoch=curriculum_epoch,
            max_curriculum_epochs=args.max_curriculum_epochs
        )
        
        train_loader = DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            shuffle=True,
            num_workers=4,
            collate_fn=collate_fn,
            pin_memory=True
        )
        
        # Train
        train_loss = train_epoch(model, train_loader, criterion, optimizer, device, epoch)
        
        # Update learning rate
        scheduler.step(train_loss)
        
        # Log
        print(f"\nEpoch {epoch+1}/{args.epochs}:")
        print(f"  Train Loss: {train_loss:.6f}")
        print(f"  Learning Rate: {optimizer.param_groups[0]['lr']:.6f}")
        print(f"  Curriculum: {curriculum_epoch}/{args.max_curriculum_epochs}")
        
        epoch_log = {
            'epoch': epoch + 1,
            'train_loss': train_loss,
            'learning_rate': optimizer.param_groups[0]['lr'],
            'curriculum_epoch': curriculum_epoch,
            'timestamp': datetime.now().isoformat()
        }
        training_log['epochs'].append(epoch_log)
        
        # Save checkpoint if best
        if train_loss < best_loss:
            best_loss = train_loss
            no_improve_epochs = 0
            
            checkpoint_path = checkpoint_dir / 'genocache_best.pt'
            model.save_checkpoint(
                checkpoint_path,
                metadata={
                    'epoch': epoch + 1,
                    'train_loss': train_loss,
                    'best_loss': best_loss,
                    'config': vars(args)
                }
            )
            
            print(f"  ✅ Saved best checkpoint (loss: {best_loss:.6f})")
        else:
            no_improve_epochs += 1
        
        # Save periodic checkpoints
        if (epoch + 1) % 5 == 0:
            periodic_path = checkpoint_dir / f'genocache_epoch{epoch+1}.pt'
            model.save_checkpoint(periodic_path, metadata=epoch_log)
            print(f"  💾 Saved epoch {epoch+1} checkpoint")
        
        # Save training log
        log_path = output_dir / 'training_log.json'
        with open(log_path, 'w') as f:
            json.dump(training_log, f, indent=2)
        
        print()
        
        # Early stopping
        if no_improve_epochs >= args.patience:
            print(f"⚠️  No improvement for {args.patience} epochs. Stopping early.")
            break
    
    # Save final checkpoint
    final_path = checkpoint_dir / 'genocache_final.pt'
    model.save_checkpoint(final_path, metadata=training_log['epochs'][-1])
    
    training_log['end_time'] = datetime.now().isoformat()
    training_log['best_loss'] = best_loss
    
    with open(log_path, 'w') as f:
        json.dump(training_log, f, indent=2)
    
    print(f"\n{'='*80}")
    print("TRAINING COMPLETE")
    print(f"{'='*80}")
    print(f"Best train loss: {best_loss:.6f}")
    print(f"Total epochs: {epoch + 1}")
    print(f"Checkpoints: {checkpoint_dir}")
    print(f"Training log: {log_path}")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
