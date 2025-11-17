#!/usr/bin/env python3
"""
GenoCache V4 - Training on chr22 (Quick Fix Phase 1)

This script trains the model on chr22 ONLY using:
1. Batch contrastive learning (InfoNCE) - NO hard negatives
2. Position shift ±51bp (exactly Lseed/10)
3. Consistent 512bp seeds throughout
4. Large batch size (8192) for diverse negatives

Timeline: 4-6 hours on 1 GPU
Expected outcome: Separation >0.6, neg_sim <0.3

Reference: NeuralAligner paper (ICLR 2026)
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from tqdm import tqdm
import time
from datetime import datetime

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder
sys.path.append(str(Path(__file__).parent))
from augmentation import DataAugmentation


class InfoNCELoss(nn.Module):
    """
    InfoNCE loss for contrastive learning
    
    Uses other samples in batch as implicit negatives.
    NO explicit hard negatives!
    
    Formula:
        L = -log(exp(sim(anchor, positive) / τ) / 
                 Σ_j exp(sim(anchor, sample_j) / τ))
    
    where sample_j includes positive and all other negatives in batch
    """
    
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
    
    def forward(self, anchor_emb, positive_emb):
        """
        Args:
            anchor_emb: [B, D] - anchor embeddings
            positive_emb: [B, D] - positive embeddings
        
        Returns:
            loss: scalar - InfoNCE loss
            pos_sim: scalar - average positive similarity
            neg_sim: scalar - average negative similarity
        """
        batch_size = anchor_emb.shape[0]
        
        # Normalize embeddings
        anchor_emb = F.normalize(anchor_emb, p=2, dim=1)
        positive_emb = F.normalize(positive_emb, p=2, dim=1)
        
        # Compute similarity matrix [B, B]
        # similarity[i,j] = cosine(anchor_i, positive_j)
        similarity = torch.mm(anchor_emb, positive_emb.T) / self.temperature
        
        # Positive similarities are on the diagonal
        pos_sim = torch.diagonal(similarity)
        
        # For each anchor, compute log-sum-exp over all positives (including diagonal)
        # This treats diagonal as positive and off-diagonal as negatives
        log_sum_exp = torch.logsumexp(similarity, dim=1)
        
        # InfoNCE loss: -log(exp(pos) / sum_exp(all))
        loss = -torch.mean(pos_sim - log_sum_exp)
        
        # Compute statistics for monitoring
        pos_sim_avg = torch.mean(pos_sim).item()
        
        # Negative similarities (off-diagonal)
        mask = torch.eye(batch_size, device=similarity.device).bool()
        neg_similarities = similarity.masked_fill(mask, float('-inf'))
        neg_sim_avg = torch.mean(neg_similarities[~mask]).item()
        
        return loss, pos_sim_avg, neg_sim_avg


def load_genome(fasta_path: Path) -> dict:
    """
    Load reference genome from FASTA file
    
    Args:
        fasta_path: Path to reference FASTA
    
    Returns:
        genome_dict: {chr_name: sequence}
    """
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
                    genome[current_chr] = ''.join(current_seq).upper()
                
                # Start new chromosome
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last chromosome
        if current_chr is not None:
            genome[current_chr] = ''.join(current_seq).upper()
    
    return genome


def generate_training_batch(genome_seq, augmenter, batch_size=8192, device='cuda'):
    """
    Generate training batch on-the-fly
    
    Args:
        genome_seq: Chromosome sequence string
        augmenter: DataAugmentation instance
        batch_size: Batch size (large for InfoNCE)
        device: torch device
    
    Returns:
        anchor_tensor: [B, 5, 512] - anchor sequences
        positive_tensor: [B, 5, 512] - positive sequences
    """
    # Generate pairs
    anchors, positives = augmenter.create_training_pairs(genome_seq, batch_size=batch_size)
    
    # Convert to tensors
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        """Convert DNA sequence to index tensor"""
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    anchor_tensors = torch.stack([seq_to_tensor(seq) for seq in anchors])
    positive_tensors = torch.stack([seq_to_tensor(seq) for seq in positives])
    
    return anchor_tensors.to(device), positive_tensors.to(device)


def compute_separation(pos_sim, neg_sim):
    """
    Compute separation metric: how well positives and negatives are separated
    
    Separation = pos_sim - neg_sim
    
    Target: >0.6 (NeuralAligner achieves >0.7)
    """
    return pos_sim - neg_sim


def train_epoch(model, genome_seq, augmenter, criterion, optimizer, device, 
                batches_per_epoch=100, batch_size=8192):
    """
    Train for one epoch
    
    Args:
        model: GenoCacheEncoder
        genome_seq: chr22 sequence
        augmenter: DataAugmentation instance
        criterion: InfoNCELoss
        optimizer: torch optimizer
        device: torch device
        batches_per_epoch: Number of batches (default 100 = ~800K examples)
        batch_size: Batch size (default 8192)
    
    Returns:
        metrics: dict with loss, pos_sim, neg_sim, separation
    """
    model.train()
    
    epoch_loss = 0
    epoch_pos_sim = 0
    epoch_neg_sim = 0
    
    pbar = tqdm(range(batches_per_epoch), desc="Training")
    for batch_idx in pbar:
        # Generate batch on-the-fly
        anchor_batch, positive_batch = generate_training_batch(
            genome_seq, augmenter, batch_size, device
        )
        
        # Forward pass
        anchor_emb = model(anchor_batch)
        positive_emb = model(positive_batch)
        
        # Compute loss
        loss, pos_sim, neg_sim = criterion(anchor_emb, positive_emb)
        
        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # Accumulate metrics
        epoch_loss += loss.item()
        epoch_pos_sim += pos_sim
        epoch_neg_sim += neg_sim
        
        # Update progress bar
        separation = compute_separation(pos_sim, neg_sim)
        pbar.set_postfix({
            'loss': f'{loss.item():.4f}',
            'sep': f'{separation:.4f}',
            'pos': f'{pos_sim:.4f}',
            'neg': f'{neg_sim:.4f}'
        })
    
    # Average metrics
    n_batches = batches_per_epoch
    return {
        'loss': epoch_loss / n_batches,
        'pos_sim': epoch_pos_sim / n_batches,
        'neg_sim': epoch_neg_sim / n_batches,
        'separation': compute_separation(epoch_pos_sim / n_batches, epoch_neg_sim / n_batches)
    }


def validate(model, genome_seq, augmenter, criterion, device, num_batches=20, batch_size=8192):
    """
    Validate on chr22
    
    Args:
        model: GenoCacheEncoder
        genome_seq: chr22 sequence
        augmenter: DataAugmentation instance
        criterion: InfoNCELoss
        device: torch device
        num_batches: Number of validation batches
        batch_size: Batch size
    
    Returns:
        metrics: dict with loss, pos_sim, neg_sim, separation
    """
    model.eval()
    
    val_loss = 0
    val_pos_sim = 0
    val_neg_sim = 0
    
    with torch.no_grad():
        for _ in range(num_batches):
            # Generate batch
            anchor_batch, positive_batch = generate_training_batch(
                genome_seq, augmenter, batch_size, device
            )
            
            # Forward pass
            anchor_emb = model(anchor_batch)
            positive_emb = model(positive_batch)
            
            # Compute loss
            loss, pos_sim, neg_sim = criterion(anchor_emb, positive_emb)
            
            # Accumulate
            val_loss += loss.item()
            val_pos_sim += pos_sim
            val_neg_sim += neg_sim
    
    # Average
    return {
        'loss': val_loss / num_batches,
        'pos_sim': val_pos_sim / num_batches,
        'neg_sim': val_neg_sim / num_batches,
        'separation': compute_separation(val_pos_sim / num_batches, val_neg_sim / num_batches)
    }


def main():
    print("=" * 80)
    print("GenoCache V4 - chr22 Training (Quick Fix)")
    print("=" * 80)
    print()
    
    # Configuration
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    CHR_NAME = "NC_000022.11"  # chr22
    CHECKPOINT_DIR = Path("/home/nebius/genocache/genocache-v4/models/checkpoints")
    LOG_DIR = Path("/home/nebius/genocache/genocache-v4/logs")
    
    BATCH_SIZE = 1024  # Large batch crucial for InfoNCE (reduced from 8192 due to GPU memory)
    BATCHES_PER_EPOCH = 200  # 200 batches × 1024 = 204K examples per epoch
    NUM_EPOCHS = 50
    LEARNING_RATE = 1e-4
    TEMPERATURE = 0.07  # InfoNCE temperature
    
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    LOG_DIR.mkdir(exist_ok=True)
    
    print("Configuration:")
    print(f"  Chromosome: {CHR_NAME}")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Batches per epoch: {BATCHES_PER_EPOCH}")
    print(f"  Examples per epoch: {BATCH_SIZE * BATCHES_PER_EPOCH:,}")
    print(f"  Epochs: {NUM_EPOCHS}")
    print(f"  Total examples: {BATCH_SIZE * BATCHES_PER_EPOCH * NUM_EPOCHS:,}")
    print(f"  Learning rate: {LEARNING_RATE}")
    print(f"  Temperature: {TEMPERATURE}")
    print()
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
        print(f"  Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print()
    
    # Load chr22
    print("Loading genome...")
    genome = load_genome(GENOME_PATH)
    
    if CHR_NAME not in genome:
        # Try alternative naming
        alt_names = ['chr22', '22', 'NC_000022.11']
        for alt in alt_names:
            if alt in genome:
                CHR_NAME = alt
                break
        else:
            print(f"❌ Could not find chr22! Available chromosomes:")
            for name in list(genome.keys())[:10]:
                print(f"   - {name}")
            return
    
    chr22_seq = genome[CHR_NAME]
    print(f"✅ Loaded {CHR_NAME}: {len(chr22_seq):,} bp")
    print()
    
    # Create augmenter
    print("Creating augmenter...")
    augmenter = DataAugmentation(seed_len=512)
    print("✅ Augmenter ready")
    print()
    
    # Create model
    print("Creating model...")
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    num_params = sum(p.numel() for p in model.parameters())
    print(f"✅ Model ready: {num_params:,} parameters")
    print()
    
    # Create loss and optimizer
    criterion = InfoNCELoss(temperature=TEMPERATURE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    
    # Training loop
    print("Starting training...")
    print()
    
    best_separation = 0
    start_time = time.time()
    
    log_file = LOG_DIR / f"training_chr22_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    with open(log_file, 'w') as f:
        f.write("epoch,loss,pos_sim,neg_sim,separation,val_loss,val_pos_sim,val_neg_sim,val_separation\n")
        
        for epoch in range(1, NUM_EPOCHS + 1):
            print(f"\n{'='*80}")
            print(f"Epoch {epoch}/{NUM_EPOCHS}")
            print(f"{'='*80}")
            
            # Train
            train_metrics = train_epoch(
                model, chr22_seq, augmenter, criterion, optimizer, device,
                batches_per_epoch=BATCHES_PER_EPOCH,
                batch_size=BATCH_SIZE
            )
            
            # Validate every 5 epochs
            if epoch % 5 == 0:
                print("\nValidating...")
                val_metrics = validate(
                    model, chr22_seq, augmenter, criterion, device,
                    num_batches=20,
                    batch_size=BATCH_SIZE
                )
                
                print(f"\nValidation Results:")
                print(f"  Loss: {val_metrics['loss']:.4f}")
                print(f"  Pos sim: {val_metrics['pos_sim']:.4f}")
                print(f"  Neg sim: {val_metrics['neg_sim']:.4f}")
                print(f"  Separation: {val_metrics['separation']:.4f}")
            else:
                val_metrics = {k: 0 for k in ['loss', 'pos_sim', 'neg_sim', 'separation']}
            
            # Print training summary
            print(f"\nTraining Summary:")
            print(f"  Loss: {train_metrics['loss']:.4f}")
            print(f"  Pos sim: {train_metrics['pos_sim']:.4f}")
            print(f"  Neg sim: {train_metrics['neg_sim']:.4f}")
            print(f"  Separation: {train_metrics['separation']:.4f}")
            
            # Log to file
            f.write(f"{epoch},{train_metrics['loss']:.6f},{train_metrics['pos_sim']:.6f},"
                   f"{train_metrics['neg_sim']:.6f},{train_metrics['separation']:.6f},"
                   f"{val_metrics['loss']:.6f},{val_metrics['pos_sim']:.6f},"
                   f"{val_metrics['neg_sim']:.6f},{val_metrics['separation']:.6f}\n")
            f.flush()
            
            # Save checkpoint if best
            if epoch % 5 == 0 and val_metrics['separation'] > best_separation:
                best_separation = val_metrics['separation']
                checkpoint_path = CHECKPOINT_DIR / f"chr22_best_sep{val_metrics['separation']:.4f}_epoch{epoch}.pt"
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'train_metrics': train_metrics,
                    'val_metrics': val_metrics,
                }, checkpoint_path)
                print(f"\n✅ Saved best checkpoint: {checkpoint_path.name}")
            
            # Also save every 10 epochs
            if epoch % 10 == 0:
                checkpoint_path = CHECKPOINT_DIR / f"chr22_epoch{epoch}.pt"
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'train_metrics': train_metrics,
                    'val_metrics': val_metrics,
                }, checkpoint_path)
                print(f"💾 Saved checkpoint: {checkpoint_path.name}")
            
            # Estimate time remaining
            elapsed = time.time() - start_time
            avg_epoch_time = elapsed / epoch
            remaining = avg_epoch_time * (NUM_EPOCHS - epoch)
            print(f"\nTime: {elapsed/3600:.1f}h elapsed, {remaining/3600:.1f}h remaining")
    
    # Final summary
    total_time = time.time() - start_time
    print("\n" + "="*80)
    print("Training Complete!")
    print("="*80)
    print(f"Total time: {total_time/3600:.2f} hours")
    print(f"Best separation: {best_separation:.4f}")
    print(f"Log file: {log_file}")
    print()
    
    # Check if we met targets
    print("Target Check:")
    if best_separation > 0.6:
        print(f"  ✅ Separation > 0.6: {best_separation:.4f}")
    else:
        print(f"  ⚠️  Separation < 0.6: {best_separation:.4f}")
    
    print()
    print("Next steps:")
    print("  1. Run validation on chr22 test reads")
    print("  2. If accuracy >80%, scale to full genome")
    print("  3. If accuracy <80%, analyze and debug")


if __name__ == "__main__":
    main()
