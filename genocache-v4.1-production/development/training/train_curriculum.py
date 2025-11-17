#!/usr/bin/env python3
"""
Curriculum Learning Training for GenoCache V4

Implements NeuralAligner's proven curriculum strategy:
- Start with clean sequences (0% errors)
- Gradually increase error rate (0% → 1% → 3% → 5% → 8% → 10%)
- Model learns robust representations that work on noisy data

This is THE KEY DIFFERENCE that makes neural alignment work!
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
from pathlib import Path
from tqdm import tqdm
import random

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder
from training.augmentation import DataAugmentation


class CurriculumGenomeDataset(Dataset):
    """
    Dataset that generates training data on-the-fly with curriculum learning
    
    Instead of pre-generating all data, we:
    1. Sample genome positions randomly
    2. Apply augmentation with CONTROLLED error rate
    3. Generate triplets in real-time
    
    This allows us to scale to 130M examples without 54GB storage!
    """
    
    def __init__(self, genome_sequence, num_examples, seed_len=512, error_rate=0.0):
        """
        Args:
            genome_sequence: Full concatenated genome string
            num_examples: Number of examples in this "epoch"
            seed_len: Seed length
            error_rate: Current curriculum error rate (0.0 to 0.10)
        """
        self.genome = genome_sequence
        self.num_examples = num_examples
        self.seed_len = seed_len
        self.error_rate = error_rate
        self.augmenter = DataAugmentation(seed_len=seed_len)
        self.base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def __len__(self):
        return self.num_examples
    
    def __getitem__(self, idx):
        """Generate one training triplet on-the-fly"""
        # Sample random position (ensure enough space)
        max_pos = len(self.genome) - self.seed_len
        anchor_pos = random.randint(0, max_pos)
        anchor_seq = self.genome[anchor_pos:anchor_pos + self.seed_len]
        
        # Skip if too many N's
        if anchor_seq.count('N') > self.seed_len * 0.1:
            # Resample
            return self.__getitem__((idx + 1) % self.num_examples)
        
        # Generate positive (augmented version with current error rate)
        positive_seq, _ = self.augmenter.augment_sequence(
            anchor_seq, 
            error_rate=self.error_rate
        )
        
        # Generate negative (random far away position)
        neg_pos = random.randint(0, max_pos)
        while abs(neg_pos - anchor_pos) < 10000:  # At least 10kb away
            neg_pos = random.randint(0, max_pos)
        negative_seq = self.genome[neg_pos:neg_pos + self.seed_len]
        
        # Convert to tensors
        anchor_tensor = self._seq_to_tensor(anchor_seq)
        positive_tensor = self._seq_to_tensor(positive_seq)
        negative_tensor = self._seq_to_tensor(negative_seq)
        
        return anchor_tensor, positive_tensor, negative_tensor
    
    def _seq_to_tensor(self, seq):
        """Convert DNA sequence to tensor (fixed length 512)"""
        if len(seq) > self.seed_len:
            seq = seq[:self.seed_len]
        elif len(seq) < self.seed_len:
            seq = seq + 'N' * (self.seed_len - len(seq))
        
        tensor = torch.zeros(self.seed_len, dtype=torch.long)
        for i, base in enumerate(seq):
            tensor[i] = self.base_to_idx.get(base.upper(), 4)
        return tensor


def train_epoch_with_error_rate(model, genome, optimizer, device, 
                                 error_rate, num_batches, batch_size):
    """
    Train one epoch with specific error rate
    
    Args:
        model: GenoCache encoder
        genome: Full genome sequence
        optimizer: Optimizer
        device: cuda/cpu
        error_rate: Error rate for this epoch (0.0-0.10)
        num_batches: Number of batches
        batch_size: Batch size
    
    Returns:
        average_loss: Average loss for epoch
    """
    model.train()
    total_loss = 0
    
    # Create dataset with current error rate
    dataset = CurriculumGenomeDataset(
        genome, 
        num_examples=num_batches * batch_size,
        error_rate=error_rate
    )
    
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,  # Already randomized in dataset
        num_workers=4,
        pin_memory=True
    )
    
    pbar = tqdm(dataloader, desc=f"Training (error_rate={error_rate:.1%})")
    for anchors, positives, negatives in pbar:
        # Move to device
        anchors = anchors.to(device)
        positives = positives.to(device)
        negatives = negatives.to(device)
        
        # Forward pass
        anchor_emb = model(anchors, use_contrastive_head=True)
        positive_emb = model(positives, use_contrastive_head=True)
        negative_emb = model(negatives, use_contrastive_head=True)
        
        # Contrastive loss
        pos_sim = F.cosine_similarity(anchor_emb, positive_emb)
        neg_sim = F.cosine_similarity(anchor_emb, negative_emb)
        loss = torch.mean(1 - pos_sim + torch.clamp(neg_sim, min=0))
        
        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        pbar.set_postfix({'loss': f'{loss.item():.4f}'})
    
    return total_loss / len(dataloader)


def main():
    print("="*80)
    print("GenoCache V4 - Curriculum Learning Training")
    print("="*80)
    
    # Configuration
    from Bio import SeqIO
    
    genome_path = '/home/nebius/genocache/GRCh38.fa'
    checkpoint_dir = Path('models/checkpoints')
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    
    batch_size = 32
    batches_per_epoch = 1000  # ~32K examples per epoch
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Curriculum schedule (NeuralAligner-style)
    curriculum = [
        {'epochs': 5, 'error_rate': 0.00, 'lr': 1e-4},  # Clean
        {'epochs': 5, 'error_rate': 0.01, 'lr': 1e-4},  # 1% error
        {'epochs': 10, 'error_rate': 0.03, 'lr': 1e-4}, # 3% error
        {'epochs': 10, 'error_rate': 0.05, 'lr': 5e-5}, # 5% error (lower LR)
        {'epochs': 15, 'error_rate': 0.08, 'lr': 5e-5}, # 8% error (ONT-like)
        {'epochs': 5, 'error_rate': 0.10, 'lr': 5e-5},  # 10% error (hard)
    ]
    
    total_epochs = sum(stage['epochs'] for stage in curriculum)
    
    print(f"\nConfiguration:")
    print(f"  Genome: {genome_path}")
    print(f"  Device: {device}")
    print(f"  Batch size: {batch_size}")
    print(f"  Batches/epoch: {batches_per_epoch}")
    print(f"  Total epochs: {total_epochs}")
    print(f"\nCurriculum:")
    for i, stage in enumerate(curriculum):
        print(f"  Stage {i+1}: {stage['epochs']} epochs @ {stage['error_rate']:.1%} error, LR={stage['lr']}")
    print()
    
    # Load genome
    print("Loading genome...")
    chromosomes = []
    for record in SeqIO.parse(genome_path, "fasta"):
        if record.id.startswith('NC_'):
            chromosomes.append(str(record.seq).upper())
            print(f"  Loaded {record.id}: {len(record.seq):,} bp")
    
    genome = 'N' * 10000.join(chromosomes)  # Concatenate with separator
    print(f"✅ Loaded genome: {len(genome):,} bp\n")
    
    # Create model
    print("Creating model...")
    model = GenoCacheEncoder(emb_dim=256, seed_len=512)
    model = model.to(device)
    print(f"✅ Model created: {model.count_parameters():,} parameters\n")
    
    # Optimizer (will adjust LR per curriculum stage)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=1e-4,
        weight_decay=1e-5
    )
    
    # Training with curriculum
    print("Starting curriculum training...\n")
    
    current_epoch = 0
    for stage_idx, stage in enumerate(curriculum):
        print("="*80)
        print(f"STAGE {stage_idx + 1}/{len(curriculum)}: "
              f"{stage['epochs']} epochs @ {stage['error_rate']:.1%} error")
        print("="*80)
        
        # Update learning rate
        for param_group in optimizer.param_groups:
            param_group['lr'] = stage['lr']
        
        # Train for this stage
        for epoch_in_stage in range(stage['epochs']):
            current_epoch += 1
            print(f"\nEpoch {current_epoch}/{total_epochs} "
                  f"(Stage {stage_idx+1}, Error={stage['error_rate']:.1%})")
            print("-"*80)
            
            # Train
            train_loss = train_epoch_with_error_rate(
                model, genome, optimizer, device,
                error_rate=stage['error_rate'],
                num_batches=batches_per_epoch,
                batch_size=batch_size
            )
            
            print(f"Loss: {train_loss:.4f}")
            
            # Save checkpoint every 5 epochs
            if current_epoch % 5 == 0:
                checkpoint_path = checkpoint_dir / f'curriculum_epoch{current_epoch}.pt'
                model.save_checkpoint(
                    checkpoint_path,
                    metadata={
                        'epoch': current_epoch,
                        'stage': stage_idx + 1,
                        'error_rate': stage['error_rate'],
                        'train_loss': train_loss
                    }
                )
                print(f"✅ Saved checkpoint")
    
    # Save final model
    final_path = checkpoint_dir / 'curriculum_final.pt'
    model.save_checkpoint(final_path, metadata={'epochs': total_epochs})
    
    print("\n" + "="*80)
    print("✅ Curriculum training complete!")
    print("="*80)


if __name__ == "__main__":
    main()
