"""
Training script for GenoCache-Align V2 - Fast CNN Tier
Multi-GPU distributed training with provenance tracking
"""

import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler
from torch.cuda.amp import autocast, GradScaler
from pathlib import Path
import numpy as np
import json
from datetime import datetime
from tqdm import tqdm
import argparse

from improved_cnn import ImprovedCNN, ContrastiveLoss, MODEL_INFO
from genomic_utils import load_training_data, SeedDataset


def setup_distributed():
    """Initialize distributed training"""
    if 'RANK' in os.environ:
        dist.init_process_group(backend='nccl')
        rank = dist.get_rank()
        world_size = dist.get_world_size()
        torch.cuda.set_device(rank)
        return rank, world_size
    else:
        return 0, 1


def train_epoch(model, dataloader, criterion, optimizer, scaler, device, epoch, rank=0):
    """Train for one epoch"""
    model.train()
    total_loss = 0
    num_batches = 0
    
    pbar = tqdm(dataloader, desc=f"Epoch {epoch}", disable=(rank != 0))
    
    for batch_idx, (seeds, positions) in enumerate(pbar):
        seeds = seeds.to(device)
        positions = positions.to(device)
        
        optimizer.zero_grad()
        
        # Mixed precision training
        with autocast():
            embeddings = model(seeds)
            loss = criterion(embeddings, positions)
        
        # Backward pass with gradient scaling
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        
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
        for seeds, positions in dataloader:
            seeds = seeds.to(device)
            positions = positions.to(device)
            
            embeddings = model(seeds)
            loss = criterion(embeddings, positions)
            
            total_loss += loss.item()
            num_batches += 1
    
    avg_loss = total_loss / num_batches
    return avg_loss


def main(args):
    # Setup distributed training
    rank, world_size = setup_distributed()
    device = torch.device(f'cuda:{rank}')
    
    if rank == 0:
        print(f"\n{'='*80}")
        print("GENOCACHE-ALIGN V2 TRAINING")
        print(f"{'='*80}")
        print(f"Model: {MODEL_INFO['name']}")
        print(f"Embedding dim: {MODEL_INFO['embedding_dim']}")
        print(f"Seed length: {MODEL_INFO['seed_length']}")
        print(f"GPUs: {world_size}")
        print(f"{'='*80}\n")
    
    # Create output directories
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    checkpoint_dir = output_dir / 'checkpoints'
    checkpoint_dir.mkdir(exist_ok=True)
    
    # Load training data
    if rank == 0:
        print("Loading training data...")
    
    seeds, positions = load_training_data(args.data_path)
    
    # Create dataset and dataloader
    dataset = SeedDataset(seeds, positions)
    
    if world_size > 1:
        sampler = DistributedSampler(
            dataset, 
            num_replicas=world_size,
            rank=rank,
            shuffle=True
        )
        shuffle = False
    else:
        sampler = None
        shuffle = True
    
    dataloader = DataLoader(dataset,
           batch_size=args.batch_size,
           shuffle=shuffle,
           sampler=sampler,
           num_workers=0,
           pin_memory=True,
           collate_fn=lambda batch: (
               torch.stack([b[0] for b in batch]),
               torch.tensor([int(b[1]) for b in batch], dtype=torch.long)
           ))
    
    if rank == 0:
        print(f"✅ Loaded {len(dataset):,} seeds")
        print(f"   Batches per epoch: {len(dataloader):,}")
        print(f"   Batch size: {args.batch_size}")
        print(f"   Effective batch size: {args.batch_size * world_size}\n")
    
    # Create model
    model = ImprovedCNN(
        out_dim=MODEL_INFO['embedding_dim'],
        input_len=MODEL_INFO['seed_length']
    ).to(device)
    
    if world_size > 1:
        model = DDP(model, device_ids=[rank])
    
    if rank == 0:
        num_params = sum(p.numel() for p in model.parameters())
        print(f"Model parameters: {num_params:,}\n")
    
    # Loss and optimizer
    criterion = ContrastiveLoss(
        temperature=args.temperature,
        margin=args.margin
    )
    
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay
    )
    
    # Learning rate scheduler
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer,
        T_max=args.epochs,
        eta_min=args.learning_rate * 0.01
    )
    
    # Gradient scaler for mixed precision
    scaler = GradScaler()
    
    # Training loop
    best_loss = float('inf')
    training_log = {
        'model': MODEL_INFO['name'],
        'config': vars(args),
        'epochs': []
    }
    
    for epoch in range(args.epochs):
        if sampler:
            sampler.set_epoch(epoch)
        
        # Train
        train_loss = train_epoch(
            model, dataloader, criterion, optimizer, 
            scaler, device, epoch, rank
        )
        
        # Update learning rate
        scheduler.step()
        
        # Log
        if rank == 0:
            print(f"\nEpoch {epoch}:")
            print(f"  Train Loss: {train_loss:.4f}")
            print(f"  Learning Rate: {scheduler.get_last_lr()[0]:.6f}")
            
            training_log['epochs'].append({
                'epoch': epoch,
                'train_loss': train_loss,
                'learning_rate': scheduler.get_last_lr()[0]
            })
            
            # Save checkpoint
            if train_loss < best_loss:
                best_loss = train_loss
                
                # Get model without DDP wrapper
                model_to_save = model.module if hasattr(model, 'module') else model
                
                checkpoint_path = checkpoint_dir / 'improved_cnn_best.pt'
                model_to_save.save_checkpoint(
                    checkpoint_path,
                    metadata={
                        'epoch': epoch,
                        'train_loss': train_loss,
                        'best_loss': best_loss,
                        'num_seeds': len(dataset),
                        'config': vars(args)
                    }
                )
                
                print(f"  ✅ Saved best checkpoint (loss: {best_loss:.4f})")
    
    # Save final checkpoint and training log
    if rank == 0:
        model_to_save = model.module if hasattr(model, 'module') else model
        
        final_path = checkpoint_dir / 'improved_cnn_final.pt'
        model_to_save.save_checkpoint(
            final_path,
            metadata={
                'epoch': args.epochs - 1,
                'train_loss': train_loss,
                'num_seeds': len(dataset),
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
        print(f"Best loss: {best_loss:.4f}")
        print(f"Checkpoints saved to: {checkpoint_dir}")
        print(f"Training log: {log_path}")
        print(f"{'='*80}\n")
    
    if world_size > 1:
        dist.destroy_process_group()


if __name__ == "__main__":
    import os
    
    parser = argparse.ArgumentParser(description='Train ImprovedCNN for genomic alignment')
    
    # Data
    parser.add_argument('--data-path', type=str, required=True,
                       help='Path to training data (.npz file)')
    parser.add_argument('--output-dir', type=str, default='genocache_v2/models/cnn_fast',
                       help='Output directory for checkpoints and logs')
    
    # Training
    parser.add_argument('--epochs', type=int, default=50,
                       help='Number of training epochs')
    parser.add_argument('--batch-size', type=int, default=256,
                       help='Batch size per GPU')
    parser.add_argument('--learning-rate', type=float, default=1e-3,
                       help='Initial learning rate')
    parser.add_argument('--weight-decay', type=float, default=1e-4,
                       help='Weight decay for regularization')
    
    # Loss
    parser.add_argument('--temperature', type=float, default=0.07,
                       help='Temperature for contrastive loss')
    parser.add_argument('--margin', type=float, default=0.5,
                       help='Margin for negative pairs')
    
    args = parser.parse_args()
    
    main(args)
