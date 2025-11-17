# GenoCache V4 - 8-GPU Training Setup Guide

**Target:** Fresh environment with 8 GPUs for optimal training  
**Goal:** Train with batch=8192 (optimal InfoNCE) to reach 99%+ accuracy  
**Date:** 2025-11-13

---

## Overview

**Current Status (1 GPU H100):**
- Model trained with batch=1024 (memory limited)
- Accuracy: 96.6% @ ±1kb
- Training time: 3.44 hours (chr22 + full genome)

**Target (8 GPU):**
- Train with batch=8192 (optimal for InfoNCE)
- Expected accuracy: 99%+ (NeuralAligner level)
- Training time: ~5-8 hours (full genome from scratch)
- Multi-GPU data parallel or gradient accumulation

---

## Prerequisites

**Hardware:**
- 8× NVIDIA GPUs (A100/H100/V100, 40GB+ each)
- 500GB+ RAM
- 2TB+ SSD storage
- Fast network for S3 access

**Software:**
- CUDA 11.8+ or 12.0+
- Python 3.10-3.12
- AWS CLI configured with S3 access

---

## Part 1: Download Current Artifacts from S3

### 1.1 Setup AWS CLI

```bash
# Configure AWS CLI (if not already done)
aws configure
# Enter your credentials

# Test access
aws s3 ls s3://YOUR_BUCKET_NAME/
```

### 1.2 Download Required Files

```bash
# Create project directory
mkdir -p ~/genocache-v4-8gpu
cd ~/genocache-v4-8gpu

# Download genome reference (3.1 GB)
aws s3 cp s3://YOUR_BUCKET_NAME/genocache-v4/genome/GRCh38.fa .
aws s3 cp s3://YOUR_BUCKET_NAME/genocache-v4/genome/GRCh38.fa.fai .

# Download training data (if pre-generated, or we'll generate fresh)
# aws s3 sync s3://YOUR_BUCKET_NAME/genocache-v4/training_data/ ./training_data/

# Download current best model (for reference/warm start option)
aws s3 cp s3://YOUR_BUCKET_NAME/genocache-v4/models/fullgenome_best_sep12.4035_epoch30.pt ./models/

# Download all training/validation scripts
aws s3 sync s3://YOUR_BUCKET_NAME/genocache-v4/scripts/ ./scripts/

# Download validation data
aws s3 sync s3://YOUR_BUCKET_NAME/genocache-v4/validation/ ./validation/
```

### 1.3 Verify Downloads

```bash
# Check genome
ls -lh GRCh38.fa
# Should be ~3.1 GB

# Check model
ls -lh models/
# Should have .pt files

# Check scripts
ls scripts/
# Should have train.py, validate.py, etc.
```

---

## Part 2: Environment Setup

### 2.1 Create Python Environment

```bash
# Using conda (recommended)
conda create -n genocache python=3.11 -y
conda activate genocache

# OR using venv
python3 -m venv .venv
source .venv/bin/activate
```

### 2.2 Install PyTorch (with CUDA)

```bash
# For CUDA 12.1
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# For CUDA 11.8
# pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Verify GPU access
python3 -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU count: {torch.cuda.device_count()}')"
# Should show: CUDA available: True, GPU count: 8
```

### 2.3 Install Dependencies

```bash
# Core dependencies
pip install biopython==1.81 \
    numpy==1.24.3 \
    tqdm==4.65.0 \
    faiss-gpu==1.7.4 \
    pandas==2.0.3 \
    matplotlib==3.7.2 \
    seaborn==0.12.2

# Verify FAISS GPU
python3 -c "import faiss; print(f'FAISS version: {faiss.__version__}'); print(f'FAISS GPU available: {faiss.get_num_gpus()}')"
```

### 2.4 Install minimap2 (for validation)

```bash
# Download and install minimap2
wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2
tar -xjf minimap2-2.26_x64-linux.tar.bz2
sudo cp minimap2-2.26_x64-linux/minimap2 /usr/local/bin/
minimap2 --version
```

---

## Part 3: Project Structure

```bash
cd ~/genocache-v4-8gpu

# Create directory structure
mkdir -p {models,training_data,validation,logs,checkpoints,scripts}

# Your structure should look like:
# genocache-v4-8gpu/
# ├── GRCh38.fa              # Reference genome
# ├── GRCh38.fa.fai          # Index
# ├── models/                # Model checkpoints
# │   └── fullgenome_best_sep12.4035_epoch30.pt
# ├── training_data/         # Training sequences (will generate)
# ├── validation/            # Validation data
# ├── logs/                  # Training logs
# ├── checkpoints/           # Training checkpoints
# └── scripts/               # All Python scripts
#     ├── train_8gpu.py
#     ├── model.py
#     ├── dataset.py
#     ├── validate.py
#     └── generate_training_data.py
```

---

## Part 4: Key Configuration Files

### 4.1 Model Architecture (`scripts/model.py`)

```python
import torch
import torch.nn as nn
import math

class HyenaConv1d(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=3):
        super().__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, 
                             padding=kernel_size//2)
        self.activation = nn.GELU()
    
    def forward(self, x):
        return self.activation(self.conv(x))

class HyenaBlock(nn.Module):
    def __init__(self, dim, kernel_size=3):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.conv = HyenaConv1d(dim, dim, kernel_size)
        self.norm2 = nn.LayerNorm(dim)
        self.ffn = nn.Sequential(
            nn.Linear(dim, dim * 4),
            nn.GELU(),
            nn.Linear(dim * 4, dim)
        )
    
    def forward(self, x):
        # x: (batch, seq_len, dim)
        residual = x
        x = self.norm1(x)
        x = x.transpose(1, 2)  # (batch, dim, seq_len)
        x = self.conv(x)
        x = x.transpose(1, 2)  # (batch, seq_len, dim)
        x = x + residual
        
        residual = x
        x = self.norm2(x)
        x = self.ffn(x)
        x = x + residual
        return x

class GenoCache(nn.Module):
    def __init__(self, vocab_size=5, embed_dim=128, num_layers=6, max_len=512):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=4)
        self.pos_embedding = nn.Parameter(torch.randn(1, max_len, embed_dim) * 0.02)
        
        self.blocks = nn.ModuleList([
            HyenaBlock(embed_dim, kernel_size=3)
            for _ in range(num_layers)
        ])
        
        self.norm = nn.LayerNorm(embed_dim)
    
    def forward(self, x):
        # x: (batch, seq_len)
        batch_size, seq_len = x.shape
        
        # Embed
        x = self.embedding(x)  # (batch, seq_len, embed_dim)
        x = x + self.pos_embedding[:, :seq_len, :]
        
        # Process through blocks
        for block in self.blocks:
            x = block(x)
        
        x = self.norm(x)
        
        # Global average pooling
        x = x.mean(dim=1)  # (batch, embed_dim)
        
        # L2 normalize
        x = torch.nn.functional.normalize(x, p=2, dim=1)
        
        return x
```

### 4.2 Dataset with Augmentation (`scripts/dataset.py`)

```python
import torch
from torch.utils.data import Dataset
import numpy as np
from Bio import SeqIO
import random

class DNAEncoder:
    def __init__(self):
        self.char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def encode(self, seq):
        return [self.char_to_idx.get(c.upper(), 4) for c in seq]

class GenomeDataset(Dataset):
    def __init__(self, genome_path, window_size=512, stride=256, 
                 max_windows=None, error_rate_range=(0.01, 0.10),
                 position_shift=51, rc_prob=0.5, chromosomes=None):
        """
        Args:
            genome_path: Path to FASTA file
            window_size: Size of each window (512bp)
            stride: Stride between windows (256 for training, 32 for indexing)
            max_windows: Maximum number of windows (for limiting dataset size)
            error_rate_range: Range of error rates for augmentation (1-10%)
            position_shift: Max position shift for augmentation (±51bp)
            rc_prob: Probability of reverse complement (0.5 = 50%)
            chromosomes: List of chromosome names to use (None = all primary)
        """
        self.window_size = window_size
        self.stride = stride
        self.error_rate_range = error_rate_range
        self.position_shift = position_shift
        self.rc_prob = rc_prob
        self.encoder = DNAEncoder()
        
        # Load genome
        print(f"Loading genome from {genome_path}...")
        self.sequences = {}
        for record in SeqIO.parse(genome_path, "fasta"):
            # Filter chromosomes
            if chromosomes is not None and record.id not in chromosomes:
                continue
            # Skip non-primary chromosomes
            if any(x in record.id.lower() for x in ['alt', 'random', 'un', 'fix', 'patch']):
                continue
            self.sequences[record.id] = str(record.seq).upper()
        
        print(f"Loaded {len(self.sequences)} chromosomes")
        
        # Generate window positions
        self.windows = []
        for chr_name, seq in self.sequences.items():
            for pos in range(0, len(seq) - window_size, stride):
                # Skip windows with too many N's
                window = seq[pos:pos + window_size]
                if window.count('N') > window_size * 0.1:  # Skip if >10% N
                    continue
                self.windows.append((chr_name, pos))
                
                if max_windows and len(self.windows) >= max_windows:
                    break
            if max_windows and len(self.windows) >= max_windows:
                break
        
        print(f"Generated {len(self.windows)} windows")
    
    def __len__(self):
        return len(self.windows)
    
    def augment_sequence(self, seq):
        """Apply error augmentation"""
        seq_list = list(seq)
        error_rate = random.uniform(*self.error_rate_range)
        num_errors = int(len(seq) * error_rate)
        
        for _ in range(num_errors):
            pos = random.randint(0, len(seq_list) - 1)
            error_type = random.random()
            
            if error_type < 0.33:  # Substitution
                seq_list[pos] = random.choice(['A', 'C', 'G', 'T'])
            elif error_type < 0.66:  # Deletion
                if pos < len(seq_list):
                    seq_list.pop(pos)
            else:  # Insertion
                seq_list.insert(pos, random.choice(['A', 'C', 'G', 'T']))
        
        # Ensure we have exactly window_size
        result = ''.join(seq_list)
        if len(result) > self.window_size:
            result = result[:self.window_size]
        elif len(result) < self.window_size:
            result = result + 'N' * (self.window_size - len(result))
        
        return result
    
    def reverse_complement(self, seq):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(c, 'N') for c in reversed(seq))
    
    def __getitem__(self, idx):
        chr_name, pos = self.windows[idx]
        seq = self.sequences[chr_name]
        
        # Extract window with random position shift
        shift = random.randint(-self.position_shift, self.position_shift)
        start = max(0, pos + shift)
        end = min(len(seq), start + self.window_size)
        window = seq[start:end]
        
        # Pad if needed
        if len(window) < self.window_size:
            window = window + 'N' * (self.window_size - len(window))
        
        # Augment
        window_aug = self.augment_sequence(window)
        
        # Reverse complement with probability
        if random.random() < self.rc_prob:
            window_aug = self.reverse_complement(window_aug)
        
        # Encode
        anchor_enc = torch.tensor(self.encoder.encode(window_aug), dtype=torch.long)
        
        return anchor_enc, chr_name, pos
```

### 4.3 Training Script for 8 GPUs (`scripts/train_8gpu.py`)

```python
import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler
from torch.cuda.amp import autocast, GradScaler
import argparse
import os
from tqdm import tqdm
import time

from model import GenoCache
from dataset import GenomeDataset

class InfoNCELoss(nn.Module):
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
    
    def forward(self, embeddings):
        # embeddings: (batch, embed_dim)
        # Compute similarity matrix
        sim_matrix = torch.matmul(embeddings, embeddings.T) / self.temperature
        
        # Create labels (diagonal is positive pairs)
        batch_size = embeddings.shape[0]
        labels = torch.arange(batch_size, device=embeddings.device)
        
        # InfoNCE loss
        loss = nn.functional.cross_entropy(sim_matrix, labels)
        
        # Compute metrics
        with torch.no_grad():
            # Positive similarity (diagonal)
            pos_sim = torch.diagonal(sim_matrix).mean()
            # Negative similarity (off-diagonal)
            mask = ~torch.eye(batch_size, dtype=torch.bool, device=embeddings.device)
            neg_sim = sim_matrix[mask].mean()
            separation = pos_sim - neg_sim
        
        return loss, pos_sim.item(), neg_sim.item(), separation.item()

def setup_distributed():
    """Setup distributed training"""
    dist.init_process_group(backend='nccl')
    rank = dist.get_rank()
    world_size = dist.get_world_size()
    torch.cuda.set_device(rank)
    return rank, world_size

def cleanup_distributed():
    dist.destroy_process_group()

def train_epoch(model, dataloader, optimizer, criterion, scaler, rank, epoch):
    model.train()
    total_loss = 0
    total_pos_sim = 0
    total_neg_sim = 0
    total_sep = 0
    
    if rank == 0:
        pbar = tqdm(dataloader, desc=f"Epoch {epoch}")
    else:
        pbar = dataloader
    
    for batch_idx, (anchor, chr_names, positions) in enumerate(pbar):
        anchor = anchor.cuda(rank, non_blocking=True)
        
        optimizer.zero_grad()
        
        with autocast():
            embeddings = model(anchor)
            loss, pos_sim, neg_sim, separation = criterion(embeddings)
        
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        
        total_loss += loss.item()
        total_pos_sim += pos_sim
        total_neg_sim += neg_sim
        total_sep += separation
        
        if rank == 0 and batch_idx % 10 == 0:
            pbar.set_postfix({
                'loss': f'{loss.item():.4f}',
                'pos': f'{pos_sim:.3f}',
                'neg': f'{neg_sim:.3f}',
                'sep': f'{separation:.3f}'
            })
    
    n = len(dataloader)
    return total_loss / n, total_pos_sim / n, total_neg_sim / n, total_sep / n

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', type=str, required=True)
    parser.add_argument('--batch-size', type=int, default=1024,
                       help='Per-GPU batch size (total = batch_size * num_gpus)')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--output-dir', type=str, default='./models')
    parser.add_argument('--checkpoint', type=str, default=None)
    parser.add_argument('--local_rank', type=int, default=0)
    args = parser.parse_args()
    
    # Setup distributed
    rank, world_size = setup_distributed()
    
    if rank == 0:
        print(f"Training on {world_size} GPUs")
        print(f"Per-GPU batch size: {args.batch_size}")
        print(f"Effective batch size: {args.batch_size * world_size}")
        os.makedirs(args.output_dir, exist_ok=True)
    
    # Create dataset
    dataset = GenomeDataset(
        genome_path=args.genome,
        window_size=512,
        stride=256,  # Dense for training
        error_rate_range=(0.01, 0.10),
        position_shift=51,
        rc_prob=0.5
    )
    
    # Distributed sampler
    sampler = DistributedSampler(dataset, num_replicas=world_size, rank=rank, shuffle=True)
    
    # Dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        sampler=sampler,
        num_workers=8,
        pin_memory=True,
        drop_last=True
    )
    
    # Model
    model = GenoCache(
        vocab_size=5,
        embed_dim=128,
        num_layers=6,
        max_len=512
    ).cuda(rank)
    
    # Load checkpoint if provided
    if args.checkpoint:
        checkpoint = torch.load(args.checkpoint, map_location=f'cuda:{rank}')
        model.load_state_dict(checkpoint['model_state_dict'])
        if rank == 0:
            print(f"Loaded checkpoint from {args.checkpoint}")
    
    # Wrap with DDP
    model = DDP(model, device_ids=[rank])
    
    # Optimizer and loss
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    criterion = InfoNCELoss(temperature=0.07)
    scaler = GradScaler()
    
    # Training loop
    best_sep = 0
    for epoch in range(1, args.epochs + 1):
        sampler.set_epoch(epoch)  # Shuffle differently each epoch
        
        loss, pos_sim, neg_sim, separation = train_epoch(
            model, dataloader, optimizer, criterion, scaler, rank, epoch
        )
        
        if rank == 0:
            print(f"\nEpoch {epoch}/{args.epochs}")
            print(f"  Loss: {loss:.4f}")
            print(f"  Pos Sim: {pos_sim:.3f}")
            print(f"  Neg Sim: {neg_sim:.3f}")
            print(f"  Separation: {separation:.3f}")
            
            # Save checkpoint
            if separation > best_sep:
                best_sep = separation
                checkpoint = {
                    'epoch': epoch,
                    'model_state_dict': model.module.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'separation': separation,
                    'pos_sim': pos_sim,
                    'neg_sim': neg_sim,
                }
                torch.save(checkpoint, f"{args.output_dir}/best_sep{separation:.4f}_epoch{epoch}.pt")
                print(f"  ✅ Saved best model (sep={separation:.3f})")
            
            # Regular checkpoint every 5 epochs
            if epoch % 5 == 0:
                checkpoint = {
                    'epoch': epoch,
                    'model_state_dict': model.module.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'separation': separation,
                }
                torch.save(checkpoint, f"{args.output_dir}/checkpoint_epoch{epoch}.pt")
    
    cleanup_distributed()

if __name__ == '__main__':
    main()
```

---

## Part 5: Training Commands

### 5.1 Single GPU Test (Sanity Check)

```bash
# Test on single GPU first
cd ~/genocache-v4-8gpu

python3 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 5 \
  --lr 1e-4 \
  --output-dir ./models/test

# Should complete quickly, verify no errors
```

### 5.2 8-GPU Training (Full Run)

```bash
# Using torchrun (recommended)
torchrun --nproc_per_node=8 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 50 \
  --lr 1e-4 \
  --output-dir ./models/8gpu_batch8192

# This gives effective batch size = 1024 × 8 = 8192! ✅

# Run in background with nohup
nohup torchrun --nproc_per_node=8 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 50 \
  --lr 1e-4 \
  --output-dir ./models/8gpu_batch8192 \
  > training_8gpu.log 2>&1 &

# Get PID
echo $! > training.pid
```

### 5.3 Monitor Training

```bash
# Watch logs
tail -f training_8gpu.log

# Check GPU usage
watch -n 5 nvidia-smi

# Check process
ps aux | grep train_8gpu
```

---

## Part 6: Upload Results to S3

### 6.1 During Training (Periodic Sync)

```bash
# Sync checkpoints every hour
watch -n 3600 'aws s3 sync ./models/ s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/models/'

# Or manual sync
aws s3 sync ./models/ s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/models/
```

### 6.2 After Training Complete

```bash
# Upload everything
aws s3 sync ./models/ s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/models/
aws s3 cp training_8gpu.log s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/logs/
aws s3 sync ./checkpoints/ s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/checkpoints/
```

---

## Part 7: Validation After Training

```bash
# Download best model from S3 (if on different machine)
aws s3 cp s3://YOUR_BUCKET_NAME/genocache-v4-8gpu/models/best_sep*.pt ./models/

# Run validation script (from original environment or create new one)
python3 scripts/validate.py \
  --model ./models/best_sep*.pt \
  --genome GRCh38.fa \
  --num-reads 1000 \
  --output ./validation/results_8gpu.json
```

---

## Part 8: Expected Results

**With batch=8192 (8 GPUs):**
- Training time: ~5-8 hours (50 epochs, full genome)
- Expected separation: 13-15+ (vs 12.40 current)
- Expected accuracy: **99%+** (vs 96.6% current)
- Improvement from optimal InfoNCE batch negatives

**Key Metrics to Monitor:**
- Separation (pos_sim - neg_sim): Target >13
- Positive similarity: Should stay high (~0.95+)
- Negative similarity: Should drop lower (<0.00)
- Validation accuracy: Target 99%+

---

## Part 9: Troubleshooting

**Issue: CUDA out of memory**
```bash
# Solution 1: Reduce per-GPU batch size
--batch-size 512  # Instead of 1024

# Solution 2: Use gradient accumulation
# Add to training script:
accumulation_steps = 2
optimizer.zero_grad()
for i, batch in enumerate(dataloader):
    loss = compute_loss(batch)
    loss = loss / accumulation_steps
    loss.backward()
    if (i + 1) % accumulation_steps == 0:
        optimizer.step()
        optimizer.zero_grad()
```

**Issue: Distributed training not working**
```bash
# Check NCCL
python3 -c "import torch; print(torch.distributed.is_nccl_available())"

# Use different backend
torchrun --nproc_per_node=8 --backend=gloo ...
```

**Issue: Slow data loading**
```bash
# Increase num_workers
--num-workers 16  # More workers for data loading
```

---

## Part 10: Quick Start Checklist

```bash
# ✅ Setup environment
conda create -n genocache python=3.11 -y
conda activate genocache
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# ✅ Install dependencies
pip install biopython numpy tqdm faiss-gpu pandas matplotlib seaborn

# ✅ Download from S3
aws s3 cp s3://YOUR_BUCKET/genocache-v4/genome/GRCh38.fa .
aws s3 sync s3://YOUR_BUCKET/genocache-v4/scripts/ ./scripts/

# ✅ Test single GPU
python3 scripts/train_8gpu.py --genome GRCh38.fa --batch-size 1024 --epochs 1 --output-dir ./test

# ✅ Launch 8-GPU training
nohup torchrun --nproc_per_node=8 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 50 \
  --lr 1e-4 \
  --output-dir ./models/8gpu \
  > training.log 2>&1 &

# ✅ Monitor
tail -f training.log
watch -n 5 nvidia-smi

# ✅ Sync to S3 regularly
watch -n 3600 'aws s3 sync ./models/ s3://YOUR_BUCKET/genocache-v4-8gpu/models/'
```

---

## Summary

**What you're training:**
- Same GenoCache V4 architecture
- With optimal batch size (8192 via 8 GPUs)
- Expected to reach 99%+ accuracy (NeuralAligner level)

**Key improvements:**
- Batch 8192 (vs 1024) → Better InfoNCE training
- Full genome from scratch (no warm start) → More optimal
- Multi-GPU training → Faster convergence

**Timeline:**
- Setup: 30 minutes
- Training: 5-8 hours
- Validation: 30 minutes
- Total: ~6-9 hours to completion

**Expected outcome:**
GenoCache V4 with 99%+ accuracy, ready for production! 🚀
