"""
GenoCache V4 - 8-GPU Distributed Training Script

Usage:
    torchrun --nproc_per_node=8 train_8gpu.py --genome GRCh38.fa --batch-size 1024 --epochs 50

This achieves effective batch size of 8192 (1024 × 8 GPUs) for optimal InfoNCE!
"""

import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler
from torch.cuda.amp import autocast, GradScaler
import argparse
import os
import time
from tqdm import tqdm
import json

from model import GenoCache
from dataset import GenomeDataset, collate_fn

class InfoNCELoss(nn.Module):
    """InfoNCE loss for contrastive learning"""
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
        self.cross_entropy = nn.CrossEntropyLoss()
    
    def forward(self, embeddings):
        """
        Args:
            embeddings: (batch, embed_dim) L2-normalized embeddings
        
        Returns:
            loss: InfoNCE loss
            pos_sim: Average positive similarity
            neg_sim: Average negative similarity
            separation: pos_sim - neg_sim
        """
        # Compute similarity matrix
        sim_matrix = torch.matmul(embeddings, embeddings.T) / self.temperature
        
        # Labels: diagonal elements are positive pairs
        batch_size = embeddings.shape[0]
        labels = torch.arange(batch_size, device=embeddings.device)
        
        # InfoNCE loss (cross-entropy)
        loss = self.cross_entropy(sim_matrix, labels)
        
        # Compute metrics (no gradient)
        with torch.no_grad():
            # Positive similarity (diagonal)
            pos_sim = torch.diagonal(sim_matrix).mean()
            
            # Negative similarity (off-diagonal)
            mask = ~torch.eye(batch_size, dtype=torch.bool, device=embeddings.device)
            neg_sim = sim_matrix[mask].mean()
            
            # Separation
            separation = pos_sim - neg_sim
        
        return loss, pos_sim.item(), neg_sim.item(), separation.item()

def setup_distributed():
    """Initialize distributed training"""
    dist.init_process_group(backend='nccl')
    rank = dist.get_rank()
    world_size = dist.get_world_size()
    torch.cuda.set_device(rank)
    return rank, world_size

def cleanup_distributed():
    """Cleanup distributed training"""
    dist.destroy_process_group()

def train_epoch(model, dataloader, optimizer, criterion, scaler, rank, epoch, args):
    """Train one epoch"""
    model.train()
    total_loss = 0
    total_pos_sim = 0
    total_neg_sim = 0
    total_sep = 0
    num_batches = 0
    
    if rank == 0:
        pbar = tqdm(dataloader, desc=f"Epoch {epoch}/{args.epochs}")
    else:
        pbar = dataloader
    
    start_time = time.time()
    
    for batch_idx, (anchor, chr_names, positions) in enumerate(pbar):
        anchor = anchor.cuda(rank, non_blocking=True)
        
        optimizer.zero_grad()
        
        # Mixed precision forward pass
        with autocast():
            embeddings = model(anchor)
            loss, pos_sim, neg_sim, separation = criterion(embeddings)
        
        # Backward pass with gradient scaling
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        
        # Accumulate metrics
        total_loss += loss.item()
        total_pos_sim += pos_sim
        total_neg_sim += neg_sim
        total_sep += separation
        num_batches += 1
        
        # Update progress bar (rank 0 only)
        if rank == 0 and batch_idx % 10 == 0:
            pbar.set_postfix({
                'loss': f'{loss.item():.4f}',
                'pos': f'{pos_sim:.3f}',
                'neg': f'{neg_sim:.3f}',
                'sep': f'{separation:.3f}'
            })
    
    epoch_time = time.time() - start_time
    
    # Average metrics
    avg_loss = total_loss / num_batches
    avg_pos_sim = total_pos_sim / num_batches
    avg_neg_sim = total_neg_sim / num_batches
    avg_sep = total_sep / num_batches
    
    return avg_loss, avg_pos_sim, avg_neg_sim, avg_sep, epoch_time

def main():
    parser = argparse.ArgumentParser(description='GenoCache V4 - 8-GPU Training')
    parser.add_argument('--genome', type=str, required=True,
                       help='Path to genome FASTA file')
    parser.add_argument('--batch-size', type=int, default=1024,
                       help='Per-GPU batch size (total = batch_size * num_gpus)')
    parser.add_argument('--epochs', type=int, default=50,
                       help='Number of training epochs')
    parser.add_argument('--lr', type=float, default=1e-4,
                       help='Learning rate')
    parser.add_argument('--weight-decay', type=float, default=0.01,
                       help='Weight decay for AdamW')
    parser.add_argument('--temperature', type=float, default=0.07,
                       help='Temperature for InfoNCE loss')
    parser.add_argument('--output-dir', type=str, default='./models',
                       help='Directory to save checkpoints')
    parser.add_argument('--checkpoint', type=str, default=None,
                       help='Path to checkpoint to resume from')
    parser.add_argument('--stride', type=int, default=256,
                       help='Stride for window extraction (256 for training)')
    parser.add_argument('--max-windows', type=int, default=None,
                       help='Maximum number of windows (None = all)')
    parser.add_argument('--num-workers', type=int, default=8,
                       help='Number of data loading workers per GPU')
    parser.add_argument('--save-every', type=int, default=5,
                       help='Save checkpoint every N epochs')
    parser.add_argument('--local_rank', type=int, default=0,
                       help='Local rank (auto-set by torchrun)')
    
    args = parser.parse_args()
    
    # Setup distributed training
    rank, world_size = setup_distributed()
    
    # Print configuration (rank 0 only)
    if rank == 0:
        print("╔" + "═" * 78 + "╗")
        print("║" + " " * 20 + "GenoCache V4 - 8-GPU Training" + " " * 29 + "║")
        print("╚" + "═" * 78 + "╝")
        print()
        print(f"Configuration:")
        print(f"  GPUs: {world_size}")
        print(f"  Per-GPU batch size: {args.batch_size}")
        print(f"  Effective batch size: {args.batch_size * world_size}")
        print(f"  Epochs: {args.epochs}")
        print(f"  Learning rate: {args.lr}")
        print(f"  Temperature: {args.temperature}")
        print(f"  Output directory: {args.output_dir}")
        print()
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Save configuration
        with open(f"{args.output_dir}/config.json", 'w') as f:
            json.dump(vars(args), f, indent=2)
    
    # Create dataset
    if rank == 0:
        print("Creating dataset...")
    
    dataset = GenomeDataset(
        genome_path=args.genome,
        window_size=512,
        stride=args.stride,
        max_windows=args.max_windows,
        error_rate_range=(0.01, 0.10),
        position_shift=51,
        rc_prob=0.5
    )
    
    # Distributed sampler
    sampler = DistributedSampler(
        dataset,
        num_replicas=world_size,
        rank=rank,
        shuffle=True,
        drop_last=True
    )
    
    # DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        sampler=sampler,
        num_workers=args.num_workers,
        pin_memory=True,
        drop_last=True,
        collate_fn=collate_fn
    )
    
    if rank == 0:
        print(f"\nDataLoader created:")
        print(f"  Total samples: {len(dataset):,}")
        print(f"  Batches per GPU: {len(dataloader):,}")
        print(f"  Samples per GPU per epoch: {len(dataloader) * args.batch_size:,}")
        print()
    
    # Create model
    model = GenoCache(
        vocab_size=5,
        embed_dim=128,
        num_layers=6,
        max_len=512
    ).cuda(rank)
    
    if rank == 0:
        print(f"Model created: {model.get_num_params():,} parameters")
    
    # Load checkpoint if provided
    start_epoch = 1
    best_sep = 0
    
    if args.checkpoint:
        checkpoint = torch.load(args.checkpoint, map_location=f'cuda:{rank}')
        model.load_state_dict(checkpoint['model_state_dict'])
        if rank == 0:
            print(f"✅ Loaded checkpoint from {args.checkpoint}")
            if 'epoch' in checkpoint:
                start_epoch = checkpoint['epoch'] + 1
                print(f"   Resuming from epoch {start_epoch}")
            if 'separation' in checkpoint:
                best_sep = checkpoint['separation']
                print(f"   Best separation so far: {best_sep:.4f}")
    
    # Wrap with DDP
    model = DDP(model, device_ids=[rank])
    
    # Optimizer and loss
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay
    )
    
    criterion = InfoNCELoss(temperature=args.temperature)
    scaler = GradScaler()
    
    if rank == 0:
        print("\nStarting training...\n")
    
    # Training history
    history = []
    
    # Training loop
    for epoch in range(start_epoch, args.epochs + 1):
        # Set epoch for sampler (important for shuffling)
        sampler.set_epoch(epoch)
        
        # Train one epoch
        loss, pos_sim, neg_sim, separation, epoch_time = train_epoch(
            model, dataloader, optimizer, criterion, scaler, rank, epoch, args
        )
        
        # Log results (rank 0 only)
        if rank == 0:
            print(f"\n{'=' * 80}")
            print(f"Epoch {epoch}/{args.epochs} - Completed in {epoch_time:.1f}s")
            print(f"{'=' * 80}")
            print(f"  Loss:       {loss:.4f}")
            print(f"  Pos Sim:    {pos_sim:.3f}")
            print(f"  Neg Sim:    {neg_sim:.3f}")
            print(f"  Separation: {separation:.3f}")
            
            # Save history
            history.append({
                'epoch': epoch,
                'loss': loss,
                'pos_sim': pos_sim,
                'neg_sim': neg_sim,
                'separation': separation,
                'time': epoch_time
            })
            
            with open(f"{args.output_dir}/history.json", 'w') as f:
                json.dump(history, f, indent=2)
            
            # Save best model
            if separation > best_sep:
                best_sep = separation
                checkpoint = {
                    'epoch': epoch,
                    'model_state_dict': model.module.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'separation': separation,
                    'pos_sim': pos_sim,
                    'neg_sim': neg_sim,
                    'loss': loss,
                    'config': vars(args)
                }
                path = f"{args.output_dir}/best_sep{separation:.4f}_epoch{epoch}.pt"
                torch.save(checkpoint, path)
                print(f"  ✅ Saved best model: {path}")
            
            # Save checkpoint periodically
            if epoch % args.save_every == 0:
                checkpoint = {
                    'epoch': epoch,
                    'model_state_dict': model.module.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'separation': separation,
                    'config': vars(args)
                }
                path = f"{args.output_dir}/checkpoint_epoch{epoch}.pt"
                torch.save(checkpoint, path)
                print(f"  💾 Saved checkpoint: {path}")
            
            print()
    
    # Final message
    if rank == 0:
        print("╔" + "═" * 78 + "╗")
        print("║" + " " * 28 + "TRAINING COMPLETE! ✅" + " " * 29 + "║")
        print("╚" + "═" * 78 + "╝")
        print()
        print(f"Best separation: {best_sep:.4f}")
        print(f"Models saved to: {args.output_dir}/")
        print()
    
    cleanup_distributed()

if __name__ == '__main__':
    main()
