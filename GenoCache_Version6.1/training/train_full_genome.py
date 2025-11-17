#!/usr/bin/env python3
"""
GenoCache v6.1 - Full Genome NAL Training

Based on successful Chr22 training (76.5% accuracy), now scaled to full genome.

Key features:
- ALL chromosomes (GRCh38)
- 24,000 batches (3x more than Chr22)
- Curriculum learning: 5% → 15% errors
- Same architecture that beat minimap2 on Chr22
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from datetime import datetime
import random

# Import encoder
sys.path.insert(0, str(Path(__file__).parent))
from encoder_nal import NALEncoder

def load_genome(fasta_path):
    """Load ALL chromosomes from GRCh38"""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous chromosome
                if current_chr is not None:
                    seq = ''.join(current_seq).upper()
                    # Filter out chromosomes with too many N's
                    if len(seq) > 1_000_000 and seq.count('N') / len(seq) < 0.1:
                        genome[current_chr] = seq
                        print(f"  Loaded {current_chr}: {len(seq):,} bp")
                
                # Start new chromosome
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last chromosome
        if current_chr is not None:
            seq = ''.join(current_seq).upper()
            if len(seq) > 1_000_000 and seq.count('N') / len(seq) < 0.1:
                genome[current_chr] = seq
                print(f"  Loaded {current_chr}: {len(seq):,} bp")
    
    total_len = sum(len(s) for s in genome.values())
    print(f"\n  Total: {len(genome)} chromosomes, {total_len:,} bp")
    return genome

def add_errors(seq, error_rate, error_types=['sub', 'ins', 'del']):
    """
    Add realistic ONT sequencing errors
    
    Implements NeuralAligner specification:
    "substitutions, insertions, and deletions applied independently 
     at each position, mimicking third-generation sequencing noise"
    """
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
                # Don't advance i (insertion)
            else:  # del
                i += 1  # Deletion (skip base)
        else:
            result.append(seq[i])
            i += 1
    
    return ''.join(result)

def generate_training_batch(genome_dict, batch_size, seed_len, error_rate, device):
    """
    Generate training batch from full genome
    
    Samples uniformly across all chromosomes based on chromosome length.
    """
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    # Weighted sampling by chromosome length
    chr_names = list(genome_dict.keys())
    chr_lengths = [len(genome_dict[c]) for c in chr_names]
    chr_weights = np.array(chr_lengths) / sum(chr_lengths)
    
    anchors = []
    positives = []
    
    for _ in range(batch_size):
        # Sample chromosome weighted by length
        chr_name = np.random.choice(chr_names, p=chr_weights)
        chr_seq = genome_dict[chr_name]
        
        # Random position
        max_start = len(chr_seq) - seed_len - 1
        if max_start < 0:
            continue
        
        start = random.randint(0, max_start)
        anchor_seq = chr_seq[start:start + seed_len]
        
        # Skip if too many N's
        if anchor_seq.count('N') > seed_len * 0.1:
            continue
        
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
    
    # Handle case where all samples were skipped
    if len(anchors) == 0:
        # Return dummy batch
        anchors = [[4] * seed_len] * batch_size
        positives = [[4] * seed_len] * batch_size
    
    return (
        torch.tensor(anchors, dtype=torch.long).to(device),
        torch.tensor(positives, dtype=torch.long).to(device)
    )

def info_nce_loss(anchor_emb, pos_emb, temperature=0.07):
    """InfoNCE contrastive loss"""
    anchor_emb = F.normalize(anchor_emb, p=2, dim=-1)
    pos_emb = F.normalize(pos_emb, p=2, dim=-1)
    
    similarity = torch.matmul(anchor_emb, pos_emb.T) / temperature
    labels = torch.arange(similarity.size(0)).to(similarity.device)
    
    loss = F.cross_entropy(similarity, labels)
    pred = similarity.argmax(dim=1)
    acc = (pred == labels).float().mean()
    
    return loss, acc

def train_full_genome(
    reference_path,
    output_dir,
    num_batches=16000,  # Conservative for GPU memory
    batch_size=32,
    emb_dim=128,
    seed_len=512,
    device='cuda'
):
    """
    Train NAL on full genome with curriculum learning
    
    Args:
        reference_path: Path to GRCh38.fa
        output_dir: Where to save models
        num_batches: 16,000 (2x more than Chr22's 8000)
        batch_size: 32 (auto-adjusts based on GPU memory)
        emb_dim: 128
        seed_len: 512bp
        device: cuda (auto-detects multi-GPU)
    """
    print("="*80)
    print("GENOCACHE v6.1 - FULL GENOME NAL TRAINING")
    print("="*80)
    
    # Auto-detect GPU resources
    n_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0
    gpu_names = []
    total_gpu_memory = 0
    
    if n_gpus > 0:
        for i in range(n_gpus):
            gpu_names.append(torch.cuda.get_device_name(i))
            total_gpu_memory += torch.cuda.get_device_properties(i).total_memory / 1e9
    
    print(f"\nGPU Resources:")
    if n_gpus > 0:
        print(f"  GPUs available: {n_gpus}")
        for i, name in enumerate(gpu_names):
            mem_gb = torch.cuda.get_device_properties(i).total_memory / 1e9
            print(f"    GPU {i}: {name} ({mem_gb:.1f} GB)")
        print(f"  Total GPU memory: {total_gpu_memory:.1f} GB")
        
        # Auto-adjust batch size based on memory
        if total_gpu_memory < 40:  # Less than 40GB total
            batch_size = max(16, batch_size // 2)
            print(f"  ⚠️  Limited GPU memory - reducing batch size to {batch_size}")
    else:
        print(f"  ⚠️  No GPU detected - training will be VERY slow!")
    
    print(f"\nConfiguration:")
    print(f"  Genome: Full GRCh38 (all chromosomes)")
    print(f"  Batches: {num_batches:,}")
    print(f"  Batch size: {batch_size}")
    print(f"  Embedding dim: {emb_dim}")
    print(f"  Seed length: {seed_len}bp")
    print(f"  Device: {device}")
    print(f"  Estimated time: 3-5 hours")
    print()
    
    # Load full genome
    genome = load_genome(reference_path)
    
    if len(genome) == 0:
        print("❌ No chromosomes loaded!")
        return None
    
    # Initialize model
    print(f"\nInitializing model...")
    model = NALEncoder(
        emb_dim=emb_dim,
        seed_len=seed_len,
        vocab_size=5,
        hidden_dim=128,
        num_layers=4
    )
    
    # Multi-GPU support
    if n_gpus > 1:
        print(f"  Using DataParallel across {n_gpus} GPUs")
        model = nn.DataParallel(model)
    
    model = model.to(device)
    
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=1e-4,
        weight_decay=0.01
    )
    
    print(f"  ✅ Model initialized ({sum(p.numel() for p in model.parameters()):,} parameters)")
    
    # Curriculum learning schedule (same as Chr22, scaled to 16K)
    print("\n" + "="*80)
    print("CURRICULUM LEARNING SCHEDULE")
    print("="*80)
    
    curriculum_schedule = [
        (0, 2000, 0.05),         # Batches 0-2000: 5% error (EASY)
        (2000, 5000, 0.08),      # Batches 2000-5000: 8% error
        (5000, 10000, 0.12),     # Batches 5000-10000: 12% error
        (10000, 16000, 0.15)     # Batches 10000+: 15% error (HARD)
    ]
    
    for start, end, rate in curriculum_schedule:
        print(f"  Batches {start:,}-{end:,}: {rate:.1%} error rate")
    
    def get_error_rate(batch_idx):
        for start, end, rate in curriculum_schedule:
            if start <= batch_idx < end:
                return rate
        return 0.15
    
    # Training loop
    print("\n" + "="*80)
    print("TRAINING")
    print("="*80)
    print()
    
    best_loss = float('inf')
    start_time = datetime.now()
    
    for batch_idx in range(num_batches):
        # Get current error rate
        error_rate = get_error_rate(batch_idx)
        
        # Generate batch
        anchor, positive = generate_training_batch(
            genome, batch_size, seed_len, error_rate, device
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
            elapsed = (datetime.now() - start_time).total_seconds()
            batches_per_sec = (batch_idx + 1) / elapsed
            remaining_batches = num_batches - (batch_idx + 1)
            eta_seconds = remaining_batches / batches_per_sec
            eta_hours = eta_seconds / 3600
            
            print(f"Batch {batch_idx+1:,}/{num_batches:,} | "
                  f"Error: {error_rate:.1%} | "
                  f"Loss: {loss.item():.4f} | "
                  f"Acc: {acc.item():.1%} | "
                  f"ETA: {eta_hours:.1f}h")
        
        # Save checkpoint every 2000 batches
        if (batch_idx + 1) % 2000 == 0:
            if loss.item() < best_loss:
                best_loss = loss.item()
                
                checkpoint_path = Path(output_dir) / f'models/genocache_nal_batch{batch_idx+1}.pt'
                torch.save({
                    'batch': batch_idx + 1,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'loss': loss.item(),
                    'acc': acc.item(),
                    'seed_len': seed_len,
                    'emb_dim': emb_dim,
                    'genome': 'GRCh38_full'
                }, checkpoint_path)
                
                print(f"  ✅ Checkpoint saved at batch {batch_idx+1:,}")
    
    # Save final model
    print("\n" + "="*80)
    print("SAVING FINAL MODEL")
    print("="*80)
    
    final_path = Path(output_dir) / f'models/genocache_nal_fullgenome_final.pt'
    torch.save({
        'batch': num_batches,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': loss.item(),
        'acc': acc.item(),
        'seed_len': seed_len,
        'emb_dim': emb_dim,
        'genome': 'GRCh38_full'
    }, final_path)
    
    print(f"  ✅ Saved final model: {final_path}")
    
    total_time = (datetime.now() - start_time).total_seconds() / 3600
    print(f"\n✅ Training complete! Total time: {total_time:.2f} hours")
    print()
    
    return model

def main():
    reference_path = '/home/nebius/genocache/GRCh38.fa'
    output_dir = '/home/nebius/genocache/GenoCache_Version6.1'
    
    # Train
    model = train_full_genome(
        reference_path=reference_path,
        output_dir=output_dir,
        num_batches=16000,  # 2x more than Chr22 (conservative for GPU)
        batch_size=32,      # Auto-adjusts based on available GPU memory
        emb_dim=128,
        seed_len=512,
        device='cuda'
    )
    
    if model:
        print("\n" + "="*80)
        print("NEXT STEPS:")
        print("="*80)
        print("1. Build full genome index: cd ../indexing && python3 build_full_genome_index.py")
        print("2. Test on synthetic reads")
        print("3. Compare with minimap2")
        print()

if __name__ == '__main__':
    main()
