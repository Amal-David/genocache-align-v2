"""
GenoCache-Align Training Pipeline
=================================

Improvements over NeurALigner training:
1. Hard negative mining from repetitive regions
2. Curriculum learning (easy -> hard)
3. Multi-task learning (contrastive + position regression)
4. Adaptive learning rate scheduling
5. Multi-GPU distributed training
6. Comprehensive logging and checkpointing
"""

import torch
import torch.nn as nn
import torch.distributed as dist
from torch.utils.data import Dataset, DataLoader
from torch.nn.parallel import DistributedDataParallel as DDP
import numpy as np
from pathlib import Path
from tqdm import tqdm
import argparse
from datetime import datetime
import json
import os

from genocache_encoder import GenoCacheEncoder, GenoCacheLoss, augment_sequence, augment_rc


class GenomicSeedDataset(Dataset):
    """
    Enhanced dataset with hard negative mining
    
    Improvements over NeurALigner's simple dataset:
    - Identifies repetitive regions for hard negative mining
    - Tracks sequence complexity scores
    - Curriculum learning support (easy samples first)
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
        self.repetitive_regions = {}  # Track repetitive regions for hard negatives
        
        print(f"\n{'='*80}")
        print("Loading Reference Genome")
        print(f"{'='*80}")
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            # Filter out sequences with too many Ns
            if seq.count('N') / len(seq) < 0.1:
                self.sequences[record.id] = seq
                # Identify repetitive regions (simple k-mer count approach)
                self.repetitive_regions[record.id] = self._find_repetitive_regions(seq)
        
        self.chrom_names = list(self.sequences.keys())
        self.chrom_lengths = [len(self.sequences[chrom]) for chrom in self.chrom_names]
        self.total_length = sum(self.chrom_lengths)
        
        print(f"✅ Loaded {len(self.sequences)} chromosomes")
        print(f"   Total length: {self.total_length:,} bp")
        print(f"   Repetitive regions: {sum(len(r) for r in self.repetitive_regions.values())} identified")
        print(f"   Hard negative ratio: {hard_negative_ratio:.2%}")
        print(f"   Curriculum epoch: {curriculum_epoch}/{max_curriculum_epochs}")
        print(f"{'='*80}\n")
    
    def _find_repetitive_regions(self, seq, kmer_size=15, threshold=10):
        """
        Identify repetitive regions by counting k-mer occurrences
        Returns list of (start, end) tuples for repetitive regions
        """
        from collections import defaultdict
        
        kmer_counts = defaultdict(int)
        kmer_positions = defaultdict(list)
        
        # Count k-mers
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i:i+kmer_size]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1
                kmer_positions[kmer].append(i)
        
        # Find repetitive k-mers
        repetitive_regions = []
        for kmer, count in kmer_counts.items():
            if count >= threshold:
                for pos in kmer_positions[kmer]:
                    # Extend region around k-mer
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
        """
        Curriculum learning: start with low error rates, gradually increase
        """
        if self.curriculum_epoch >= self.max_curriculum_epochs:
            # Full difficulty
            return np.random.uniform(0.01, 0.10)
        else:
            # Gradually increase difficulty
            progress = self.curriculum_epoch / self.max_curriculum_epochs
            min_error = 0.01
            max_error = 0.01 + (0.10 - 0.01) * progress
            return np.random.uniform(min_error, max_error)
    
    def _sample_hard_negative(self, chrom, exclude_start, exclude_end):
        """
        Sample a hard negative from repetitive regions
        Avoid sampling near the anchor position
        """
        repetitive = self.repetitive_regions[chrom]
        
        if not repetitive:
            # No repetitive regions, sample randomly
            return self._sample_random_seed(chrom)
        
        # Sample from repetitive region
        region = repetitive[np.random.randint(0, len(repetitive))]
        region_start, region_end = region
        
        # Ensure not too close to anchor
        if region_end - region_start > self.seed_len * 3:
            # Sample within region, avoiding anchor
            if region_start < exclude_start - self.seed_len * 2:
                pos = np.random.randint(region_start, exclude_start - self.seed_len * 2)
            elif region_end > exclude_end + self.seed_len * 2:
                pos = np.random.randint(exclude_end + self.seed_len * 2, 
                                       region_end - self.seed_len)
            else:
                pos = np.random.randint(region_start, region_end - self.seed_len)
        else:
            # Region too small, just sample from it
            pos = np.random.randint(region_start, max(region_start + 1, region_end - self.seed_len))
        
        seq = self.sequences[chrom][pos:pos + self.seed_len]
        
        # Skip if too many Ns
        if seq.count('N') > self.seed_len * 0.1:
            return self._sample_random_seed(chrom)
        
        return seq
    
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
            return self.__getitem__(idx)  # Retry with different chromosome
        
        pos = np.random.randint(0, max_pos)
        seed = chrom_seq[pos:pos + self.seed_len]
        
        # Skip if too many Ns
        if seed.count('N') > self.seed_len * 0.1:
            return self.__getitem__(idx)
        
        # Generate augmented version (positive pair)
        error_rate = self._get_curriculum_error_rate()
        augmented = augment_sequence(seed, error_rate=error_rate, seed_len=self.seed_len)
        
        # RC augmentation (apply to both with 50% probability each)
        seed = augment_rc(seed, prob=0.5)
        augmented = augment_rc(augmented, prob=0.5)
        
        # Sample hard negatives
        hard_negatives = []
        if np.random.random() < self.hard_negative_ratio:
            # Sample 2-4 hard negatives
            num_hard_negs = np.random.randint(2, 5)
            for _ in range(num_hard_negs):
                hard_neg = self._sample_hard_negative(chrom, pos, pos + self.seed_len)
                if hard_neg is not None:
                    hard_neg = augment_rc(hard_neg, prob=0.5)  # RC augment hard negatives too
                    hard_negatives.append(hard_neg)
        
        return seed, augmented, hard_negatives


def collate_fn(batch):
    """
    Convert DNA strings to tensors
    Handle variable number of hard negatives
    """
    seeds, augmented, hard_negatives_list = zip(*batch)
    
    # Convert to tensors
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def str_to_tensor(seq):
        return torch.tensor([mapping.get(base, 4) for base in seq.upper()], dtype=torch.long)
    
    seed_tensors = torch.stack([str_to_tensor(s) for s in seeds])
    aug_tensors = torch.stack([str_to_tensor(a) for a in augmented])
    
    # Handle hard negatives (variable number per sample)
    max_hard_negs = max(len(hn) for hn in hard_negatives_list)
    if max_hard_negs > 0:
        # Pad to max number of hard negatives
        hard_neg_tensors = []
        for hard_negs in hard_negatives_list:
            if hard_negs:
                tensors = [str_to_tensor(hn) for hn in hard_negs]
                # Pad if needed
                while len(tensors) < max_hard_negs:
                    tensors.append(tensors[0].clone())  # Duplicate first one as padding
                hard_neg_tensors.append(torch.stack(tensors))
            else:
                # No hard negatives, use random
                hard_neg_tensors.append(seed_tensors[:max_hard_negs].clone())
        hard_neg_tensors = torch.stack(hard_neg_tensors)
    else:
        hard_neg_tensors = None
    
    return seed_tensors, aug_tensors, hard_neg_tensors


def train_epoch(model, dataloader, criterion, optimizer, device, epoch, rank=0):
    """Train for one epoch with distributed training support"""
    model.train()
    total_loss = 0
    num_batches = 0
    
    if rank == 0:
        pbar = tqdm(dataloader, desc=f"Epoch {epoch}")
    else:
        pbar = dataloader
    
    for batch_idx, (seeds, augmented, hard_negatives) in enumerate(pbar):
        seeds = seeds.to(device)
        augmented = augmented.to(device)
        if hard_negatives is not None:
            hard_negatives = hard_negatives.to(device)
        
        optimizer.zero_grad()
        
        # Forward pass
        embeddings = model(seeds, use_contrastive_head=True)
        aug_embeddings = model(augmented, use_contrastive_head=True)
        
        # Encode hard negatives if present
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
        
        if rank == 0:
            pbar.set_postfix({'loss': f'{loss.item():.4f}'})
    
    avg_loss = total_loss / num_batches
    return avg_loss


def validate(model, dataloader, criterion, device, rank=0):
    """Validate model"""
    model.eval()
    total_loss = 0
    num_batches = 0
    
    with torch.no_grad():
        for seeds, augmented, hard_negatives in dataloader:
            seeds = seeds.to(device)
            augmented = augmented.to(device)
            if hard_negatives is not None:
                hard_negatives = hard_negatives.to(device)
            
            embeddings = model(seeds, use_contrastive_head=True)
            aug_embeddings = model(augmented, use_contrastive_head=True)
            
            # Encode hard negatives
            hard_neg_embeddings = None
            if hard_negatives is not None:
                B, K, L = hard_negatives.shape
                hard_negatives_flat = hard_negatives.view(B * K, L)
                hard_neg_embeddings_flat = model(hard_negatives_flat, use_contrastive_head=True)
                hard_neg_embeddings = hard_neg_embeddings_flat.view(B, K, -1)
            
            loss = criterion(embeddings, aug_embeddings, hard_neg_embeddings)
            
            total_loss += loss.item()
            num_batches += 1
    
    avg_loss = total_loss / num_batches
    return avg_loss


def setup_distributed(rank, world_size):
    """Setup distributed training"""
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12355'
    dist.init_process_group("nccl", rank=rank, world_size=world_size)


def cleanup_distributed():
    """Cleanup distributed training"""
    dist.destroy_process_group()


def main_worker(rank, world_size, args):
    """Main training function for each process"""
    
    if world_size > 1:
        setup_distributed(rank, world_size)
    
    device = torch.device(f'cuda:{rank}')
    
    if rank == 0:
        print(f"\n{'='*80}")
        print("GENOCACHE-ALIGN ENCODER TRAINING")
        print(f"{'='*80}")
        print(f"Seed length: {args.seed_len}")
        print(f"Embedding dim: {args.emb_dim}")
        print(f"Devices: {world_size} GPUs")
        print(f"Batch size per GPU: {args.batch_size}")
        print(f"Effective batch size: {args.batch_size * world_size}")
        print(f"Hard negative mining: {args.hard_negative_ratio:.2%}")
        print(f"Curriculum learning: {args.max_curriculum_epochs} epochs")
        print(f"{'='*80}\n")
    
    # Create output directories
    if rank == 0:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        checkpoint_dir = output_dir / 'checkpoints'
        checkpoint_dir.mkdir(exist_ok=True)
    
    # Create model
    model = GenoCacheEncoder(
        emb_dim=args.emb_dim,
        seed_len=args.seed_len
    ).to(device)
    
    if world_size > 1:
        model = DDP(model, device_ids=[rank])
    
    if rank == 0:
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
    
    # Learning rate scheduler
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=3,
        verbose=(rank == 0)
    )
    
    # Training loop
    best_loss = float('inf')
    no_improve_epochs = 0
    
    training_log = {
        'model': 'GenoCacheEncoder',
        'config': vars(args),
        'epochs': []
    }
    
    for epoch in range(args.epochs):
        # Update curriculum
        curriculum_epoch = min(epoch, args.max_curriculum_epochs)
        
        # Create datasets for this epoch
        train_dataset = GenomicSeedDataset(
            args.fasta,
            seed_len=args.seed_len,
            samples_per_epoch=args.samples_per_epoch,
            hard_negative_ratio=args.hard_negative_ratio,
            curriculum_epoch=curriculum_epoch,
            max_curriculum_epochs=args.max_curriculum_epochs
        )
        
        val_dataset = GenomicSeedDataset(
            args.fasta,
            seed_len=args.seed_len,
            samples_per_epoch=args.samples_per_epoch // 10,
            hard_negative_ratio=args.hard_negative_ratio,
            curriculum_epoch=curriculum_epoch,
            max_curriculum_epochs=args.max_curriculum_epochs
        )
        
        # Create data loaders
        train_sampler = torch.utils.data.distributed.DistributedSampler(
            train_dataset, num_replicas=world_size, rank=rank
        ) if world_size > 1 else None
        
        train_loader = DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            shuffle=(train_sampler is None),
            sampler=train_sampler,
            num_workers=4,
            collate_fn=collate_fn,
            pin_memory=True
        )
        
        val_loader = DataLoader(
            val_dataset,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=4,
            collate_fn=collate_fn,
            pin_memory=True
        )
        
        # Train
        train_loss = train_epoch(
            model, train_loader, criterion, optimizer, device, epoch, rank
        )
        
        # Validate
        val_loss = validate(model, val_loader, criterion, device, rank)
        
        # Synchronize validation loss across GPUs
        if world_size > 1:
            val_loss_tensor = torch.tensor(val_loss).to(device)
            dist.all_reduce(val_loss_tensor)
            val_loss = val_loss_tensor.item() / world_size
        
        # Update learning rate
        scheduler.step(val_loss)
        
        # Log
        if rank == 0:
            print(f"\nEpoch {epoch}:")
            print(f"  Train Loss: {train_loss:.4f}")
            print(f"  Val Loss: {val_loss:.4f}")
            print(f"  Learning Rate: {optimizer.param_groups[0]['lr']:.6f}")
            print(f"  Curriculum: {curriculum_epoch}/{args.max_curriculum_epochs}")
            
            training_log['epochs'].append({
                'epoch': epoch,
                'train_loss': train_loss,
                'val_loss': val_loss,
                'learning_rate': optimizer.param_groups[0]['lr'],
                'curriculum_epoch': curriculum_epoch
            })
            
            # Save checkpoint if best
            if val_loss < best_loss:
                best_loss = val_loss
                no_improve_epochs = 0
                
                checkpoint_path = checkpoint_dir / 'genocache_best.pt'
                if world_size > 1:
                    model.module.save_checkpoint(
                        checkpoint_path,
                        metadata={
                            'epoch': epoch,
                            'train_loss': train_loss,
                            'val_loss': val_loss,
                            'best_loss': best_loss,
                            'config': vars(args)
                        }
                    )
                else:
                    model.save_checkpoint(
                        checkpoint_path,
                        metadata={
                            'epoch': epoch,
                            'train_loss': train_loss,
                            'val_loss': val_loss,
                            'best_loss': best_loss,
                            'config': vars(args)
                        }
                    )
                
                print(f"  ✅ Saved best checkpoint (val_loss: {best_loss:.4f})")
            else:
                no_improve_epochs += 1
            
            # Early stopping
            if no_improve_epochs >= args.patience:
                print(f"\n⚠️  No improvement for {args.patience} epochs. Stopping early.")
                break
    
    # Save final checkpoint
    if rank == 0:
        final_path = checkpoint_dir / 'genocache_final.pt'
        if world_size > 1:
            model.module.save_checkpoint(
                final_path,
                metadata={
                    'epoch': epoch,
                    'train_loss': train_loss,
                    'val_loss': val_loss,
                    'config': vars(args)
                }
            )
        else:
            model.save_checkpoint(
                final_path,
                metadata={
                    'epoch': epoch,
                    'train_loss': train_loss,
                    'val_loss': val_loss,
                    'config': vars(args)
                }
            )
        
        # Save training log
        log_path = output_dir / 'training_log.json'
        with open(log_path, 'w') as f:
            json.dump(training_log, f, indent=2)
        
        print(f"\n{'='*80}")
        print("TRAINING COMPLETE")
        print(f"{'='*80}")
        print(f"Best val loss: {best_loss:.4f}")
        print(f"Checkpoints: {checkpoint_dir}")
        print(f"Training log: {log_path}")
        print(f"{'='*80}\n")
    
    if world_size > 1:
        cleanup_distributed()


def main(args):
    """Main entry point"""
    world_size = torch.cuda.device_count()
    
    if world_size > 1:
        print(f"🚀 Launching distributed training on {world_size} GPUs")
        torch.multiprocessing.spawn(
            main_worker,
            args=(world_size, args),
            nprocs=world_size,
            join=True
        )
    else:
        print(f"🚀 Launching single-GPU training")
        main_worker(0, 1, args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train GenoCacheEncoder')
    
    # Data
    parser.add_argument('--fasta', type=str, required=True,
                       help='Path to reference genome FASTA')
    parser.add_argument('--output-dir', type=str, default='genocache_models',
                       help='Output directory')
    
    # Model
    parser.add_argument('--seed-len', type=int, default=512,
                       help='Seed length')
    parser.add_argument('--emb-dim', type=int, default=256,
                       help='Embedding dimension (2× NeurALigner)')
    
    # Training
    parser.add_argument('--epochs', type=int, default=50,
                       help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=2048,
                       help='Batch size per GPU')
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
    
    # Advanced features
    parser.add_argument('--hard-negative-ratio', type=float, default=0.3,
                       help='Ratio of samples with hard negatives')
    parser.add_argument('--hard-negative-weight', type=float, default=2.0,
                       help='Weight for hard negative loss')
    parser.add_argument('--max-curriculum-epochs', type=int, default=10,
                       help='Number of curriculum learning epochs')
    
    args = parser.parse_args()
    
    main(args)