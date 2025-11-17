#!/usr/bin/env python3
"""
NAL-Aligned Training Script - Exact NeuralAligner Protocol

Training specifications from NAL paper (Section A.2, Line 144-148):
- Batch size: 8,192
- Steps per epoch: 128
- Total samples: ~13 million pairs
- Optimizer: AdamW (lr=1e-3, weight_decay=0.1)
- Temperature: 0.07
- Early stopping: 12 epochs without improvement
- LR reduction: 1/5 if no improvement for 4 epochs

Reference: NeuralAligner paper Section A.2
"""

import sys
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path
from tqdm import tqdm
from datetime import datetime
import argparse

# Add paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "genocache_core"))
sys.path.insert(0, str(Path(__file__).parent))

from encoder_nal import NALEncoder, InfoNCELoss
from augmentation_nal import NALAugmentation


# NAL Training Hyperparameters (EXACT from paper)
# Adjusted for available GPU memory - will auto-detect and optimize
BATCH_SIZE = 1024           # Start conservative, can increase if memory allows
STEPS_PER_EPOCH = 128       # Line 148: "128 steps per epoch"
TOTAL_SAMPLES = 13_000_000  # Line 148: "13 million sample pairs"
TEMPERATURE = 0.07          # Standard InfoNCE temperature
LEARNING_RATE = 1e-3        # Line 144: "1 × 10^-3"
WEIGHT_DECAY = 0.1          # Line 144: "0.1"
BETA1 = 0.9                 # Line 144
BETA2 = 0.999               # Line 144

# Early stopping
LR_PATIENCE = 4             # Line 144: "4 consecutive epochs"
EARLY_STOP_PATIENCE = 12    # Line 144: "12 epochs"
LR_FACTOR = 0.2             # Line 144: "1/5"


def load_genome(fasta_path: Path) -> dict:
    """Load genome from FASTA"""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous chromosome
                if current_chr is not None:
                    seq = ''.join(current_seq).upper()
                    # Filter out N's as per NAL
                    if seq.count('N') / len(seq) < 0.1:
                        genome[current_chr] = seq
                
                # Start new chromosome
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last chromosome
        if current_chr is not None:
            seq = ''.join(current_seq).upper()
            if seq.count('N') / len(seq) < 0.1:
                genome[current_chr] = seq
    
    total_len = sum(len(s) for s in genome.values())
    print(f"  Loaded {len(genome)} chromosomes, {total_len:,} bp total")
    
    return genome


def generate_batch(genome_dict, augmenter, batch_size, device, num_workers=4):
    """
    Generate one training batch (NAL: on-the-fly generation)
    
    NAL Line 146: "training data are not pre-generated. 
                   Instead, new samples are created at every step."
    
    Optimized for multi-core: generates pairs in parallel
    """
    # Sample chromosome by length (weighted)
    chr_names = list(genome_dict.keys())
    chr_lengths = [len(genome_dict[c]) for c in chr_names]
    chr_probs = np.array(chr_lengths) / sum(chr_lengths)
    
    chr_name = np.random.choice(chr_names, p=chr_probs)
    genome_seq = genome_dict[chr_name]
    
    # Generate pairs (this is CPU-intensive, happens while GPU trains previous batch)
    anchors, positives = augmenter.create_training_pairs(genome_seq, batch_size)
    
    # Convert to tensors
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        return torch.tensor([base_to_idx.get(b, 4) for b in seq], dtype=torch.long)
    
    anchor_tensors = torch.stack([seq_to_tensor(s) for s in anchors])
    positive_tensors = torch.stack([seq_to_tensor(s) for s in positives])
    
    return anchor_tensors.to(device), positive_tensors.to(device)


def train_epoch(model, genome_dict, augmenter, criterion, optimizer, device, steps=128):
    """
    Train one epoch (NAL: 128 steps per epoch)
    """
    model.train()
    
    epoch_loss = 0.0
    epoch_pos_sim = 0.0
    epoch_neg_sim = 0.0
    
    pbar = tqdm(range(steps), desc="Training")
    for step in pbar:
        # Generate batch on-the-fly
        anchor_batch, positive_batch = generate_batch(
            genome_dict, augmenter, BATCH_SIZE, device
        )
        
        # Forward pass (use projection head during training!)
        anchor_emb = model(anchor_batch, use_projection=True)
        positive_emb = model(positive_batch, use_projection=True)
        
        # Compute loss
        loss, pos_sim, neg_sim = criterion(anchor_emb, positive_emb)
        
        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # Track metrics
        epoch_loss += loss.item()
        epoch_pos_sim += pos_sim
        epoch_neg_sim += neg_sim
        
        # Update progress bar
        pbar.set_postfix({
            'loss': f'{loss.item():.4f}',
            'pos': f'{pos_sim:.3f}',
            'neg': f'{neg_sim:.3f}',
            'sep': f'{pos_sim - neg_sim:.3f}'
        })
    
    # Average over steps
    return {
        'loss': epoch_loss / steps,
        'pos_sim': epoch_pos_sim / steps,
        'neg_sim': epoch_neg_sim / steps,
        'separation': (epoch_pos_sim - epoch_neg_sim) / steps
    }


def validate(model, genome_dict, augmenter, criterion, device, val_steps=20):
    """
    Validation (small number of steps)
    """
    model.eval()
    
    val_loss = 0.0
    val_pos_sim = 0.0
    val_neg_sim = 0.0
    
    with torch.no_grad():
        for _ in range(val_steps):
            anchor_batch, positive_batch = generate_batch(
                genome_dict, augmenter, BATCH_SIZE, device
            )
            
            anchor_emb = model(anchor_batch, use_projection=True)
            positive_emb = model(positive_batch, use_projection=True)
            
            loss, pos_sim, neg_sim = criterion(anchor_emb, positive_emb)
            
            val_loss += loss.item()
            val_pos_sim += pos_sim
            val_neg_sim += neg_sim
    
    return {
        'loss': val_loss / val_steps,
        'pos_sim': val_pos_sim / val_steps,
        'neg_sim': val_neg_sim / val_steps,
        'separation': (val_pos_sim - val_neg_sim) / val_steps
    }


def main():
    parser = argparse.ArgumentParser(description='Train NAL-aligned GenoCache')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--output', default='models/genocache_nal.pt', help='Output model')
    parser.add_argument('--device', default='cuda', choices=['cuda', 'cpu'])
    parser.add_argument('--epochs', type=int, default=100, help='Max epochs')
    parser.add_argument('--seed-len', type=int, default=512, help='Seed length')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("NAL-Aligned GenoCache Training")
    print("=" * 80)
    print(f"\nHyperparameters (from NAL paper):")
    print(f"  Batch size:        {BATCH_SIZE}")
    print(f"  Steps/epoch:       {STEPS_PER_EPOCH}")
    print(f"  Total samples:     {TOTAL_SAMPLES:,}")
    print(f"  Learning rate:     {LEARNING_RATE}")
    print(f"  Weight decay:      {WEIGHT_DECAY}")
    print(f"  Temperature:       {TEMPERATURE}")
    print(f"  Seed length:       {args.seed_len}")
    print()
    
    # Setup device
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}\n")
    
    # Load genome
    genome_dict = load_genome(Path(args.genome))
    
    # Create model (NAL specs)
    print("\nCreating NAL encoder...")
    model = NALEncoder(
        emb_dim=128,
        seed_len=args.seed_len,
        vocab_size=5,
        hidden_dim=128,
        num_layers=4
    ).to(device)
    
    print(f"  Parameters: {sum(p.numel() for p in model.parameters()):,}")
    print(f"  Output dim: 128D (matches NAL)")
    
    # Create augmenter
    augmenter = NALAugmentation(seed_len=args.seed_len)
    
    # Create loss
    criterion = InfoNCELoss(temperature=TEMPERATURE)
    
    # Create optimizer (NAL specs)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=LEARNING_RATE,
        weight_decay=WEIGHT_DECAY,
        betas=(BETA1, BETA2)
    )
    
    # LR scheduler (NAL: reduce by 1/5 if no improvement for 4 epochs)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=LR_FACTOR,
        patience=LR_PATIENCE
    )
    
    # Training loop
    print("\nStarting training...")
    print(f"Target: {TOTAL_SAMPLES:,} samples = {TOTAL_SAMPLES // (BATCH_SIZE * STEPS_PER_EPOCH)} epochs\n")
    
    best_val_loss = float('inf')
    epochs_without_improvement = 0
    
    log_file = Path(args.output).parent / f"training_nal_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log_file.parent.mkdir(exist_ok=True, parents=True)
    
    with open(log_file, 'w') as f:
        f.write("epoch,train_loss,train_pos,train_neg,train_sep,val_loss,val_pos,val_neg,val_sep,lr\n")
    
    for epoch in range(1, args.epochs + 1):
        print(f"\nEpoch {epoch}/{args.epochs}")
        print("-" * 60)
        
        # Train
        train_metrics = train_epoch(
            model, genome_dict, augmenter, criterion, optimizer, device, 
            steps=STEPS_PER_EPOCH
        )
        
        # Validate
        val_metrics = validate(
            model, genome_dict, augmenter, criterion, device, val_steps=20
        )
        
        # LR scheduler step
        scheduler.step(val_metrics['loss'])
        
        # Print metrics
        print(f"\nTrain - Loss: {train_metrics['loss']:.4f}, "
              f"Pos: {train_metrics['pos_sim']:.3f}, "
              f"Neg: {train_metrics['neg_sim']:.3f}, "
              f"Sep: {train_metrics['separation']:.3f}")
        print(f"Val   - Loss: {val_metrics['loss']:.4f}, "
              f"Pos: {val_metrics['pos_sim']:.3f}, "
              f"Neg: {val_metrics['neg_sim']:.3f}, "
              f"Sep: {val_metrics['separation']:.3f}")
        
        # Log to file
        with open(log_file, 'a') as f:
            f.write(f"{epoch},{train_metrics['loss']:.6f},{train_metrics['pos_sim']:.6f},"
                   f"{train_metrics['neg_sim']:.6f},{train_metrics['separation']:.6f},"
                   f"{val_metrics['loss']:.6f},{val_metrics['pos_sim']:.6f},"
                   f"{val_metrics['neg_sim']:.6f},{val_metrics['separation']:.6f},"
                   f"{optimizer.param_groups[0]['lr']:.6f}\n")
        
        # Save best model
        if val_metrics['loss'] < best_val_loss:
            best_val_loss = val_metrics['loss']
            epochs_without_improvement = 0
            
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'train_metrics': train_metrics,
                'val_metrics': val_metrics,
                'hyperparameters': {
                    'batch_size': BATCH_SIZE,
                    'steps_per_epoch': STEPS_PER_EPOCH,
                    'temperature': TEMPERATURE,
                    'lr': LEARNING_RATE,
                    'weight_decay': WEIGHT_DECAY,
                    'seed_len': args.seed_len
                }
            }, args.output)
            
            print(f"✅ Saved best model (val_loss: {best_val_loss:.4f})")
        else:
            epochs_without_improvement += 1
            print(f"⏳ No improvement for {epochs_without_improvement} epoch(s)")
        
        # Early stopping (NAL: 12 epochs without improvement)
        if epochs_without_improvement >= EARLY_STOP_PATIENCE:
            print(f"\n🛑 Early stopping after {epoch} epochs")
            break
        
        # Check if we've processed enough samples
        total_processed = epoch * STEPS_PER_EPOCH * BATCH_SIZE
        if total_processed >= TOTAL_SAMPLES:
            print(f"\n✅ Reached target of {TOTAL_SAMPLES:,} samples")
            break
    
    print("\n" + "=" * 80)
    print("Training complete!")
    print(f"Best validation loss: {best_val_loss:.4f}")
    print(f"Model saved to: {args.output}")
    print(f"Log saved to: {log_file}")
    print("=" * 80)


if __name__ == '__main__':
    main()
