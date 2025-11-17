#!/usr/bin/env python3
"""
Training on 10M examples - improved version
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import h5py
import numpy as np
from pathlib import Path
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder


class HDF5Dataset(Dataset):
    """Dataset loader for HDF5 training data"""
    
    def __init__(self, h5_path):
        self.h5_path = h5_path
        with h5py.File(h5_path, 'r') as f:
            self.length = len(f['anchors'])
        
        self.base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, idx):
        with h5py.File(self.h5_path, 'r') as f:
            anchor = f['anchors'][idx].decode('utf-8')
            positive = f['positives'][idx].decode('utf-8')
            negative = f['negatives'][idx].decode('utf-8')
        
        anchor_tensor = self._seq_to_tensor(anchor)
        positive_tensor = self._seq_to_tensor(positive)
        negative_tensor = self._seq_to_tensor(negative)
        
        return anchor_tensor, positive_tensor, negative_tensor
    
    def _seq_to_tensor(self, seq):
        """Convert DNA sequence to tensor (fixed length 512)"""
        if len(seq) > 512:
            seq = seq[:512]
        elif len(seq) < 512:
            seq = seq + 'N' * (512 - len(seq))
        
        tensor = torch.zeros(512, dtype=torch.long)
        for i, base in enumerate(seq):
            tensor[i] = self.base_to_idx.get(base.upper(), 4)
        return tensor


def train_epoch(model, dataloader, optimizer, device):
    """Train for one epoch"""
    model.train()
    total_loss = 0
    
    pbar = tqdm(dataloader, desc="Training")
    for batch_idx, (anchors, positives, negatives) in enumerate(pbar):
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
        pbar.set_postfix({'loss': loss.item()})
    
    return total_loss / len(dataloader)


def validate(model, dataloader, device):
    """Validation pass"""
    model.eval()
    pos_sims = []
    neg_sims = []
    
    with torch.no_grad():
        for anchors, positives, negatives in tqdm(dataloader, desc="Validating"):
            anchors = anchors.to(device)
            positives = positives.to(device)
            negatives = negatives.to(device)
            
            anchor_emb = model(anchors)
            positive_emb = model(positives)
            negative_emb = model(negatives)
            
            pos_sim = F.cosine_similarity(anchor_emb, positive_emb)
            neg_sim = F.cosine_similarity(anchor_emb, negative_emb)
            
            pos_sims.extend(pos_sim.cpu().numpy())
            neg_sims.extend(neg_sim.cpu().numpy())
    
    return np.mean(pos_sims), np.mean(neg_sims)


def main():
    print("="*80)
    print("GenoCache V4 - Training on 10M Examples")
    print("="*80)
    
    # Configuration
    data_path = 'data/training_10M.h5'
    checkpoint_dir = Path('models/checkpoints')
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    
    batch_size = 64  # Increased from 32
    num_epochs = 10  # More epochs for larger dataset
    learning_rate = 1e-4
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\nConfiguration:")
    print(f"  Data: {data_path}")
    print(f"  Device: {device}")
    print(f"  Batch size: {batch_size}")
    print(f"  Epochs: {num_epochs}")
    print(f"  Learning rate: {learning_rate}")
    print()
    
    # Load dataset
    print("Loading dataset...")
    dataset = HDF5Dataset(data_path)
    train_size = int(0.9 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(
        dataset, [train_size, val_size]
    )
    
    train_loader = DataLoader(
        train_dataset, 
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        pin_memory=True
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=2
    )
    
    print(f"✅ Loaded dataset: {len(train_dataset):,} train, {len(val_dataset):,} val\n")
    
    # Create model (or load from 1M checkpoint)
    print("Creating model...")
    model = GenoCacheEncoder(emb_dim=256, seed_len=512)
    
    # Try to load from 1M checkpoint as starting point
    checkpoint_1M = Path('models/checkpoints/best_model.pt')
    if checkpoint_1M.exists():
        print("  Loading weights from 1M training...")
        checkpoint = torch.load(checkpoint_1M, map_location=device)
        model.load_state_dict(checkpoint['model_state_dict'])
        print("  ✅ Loaded pretrained weights (warm start)")
    
    model = model.to(device)
    print(f"✅ Model ready: {model.count_parameters():,} parameters\n")
    
    # Optimizer
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        weight_decay=1e-5
    )
    
    # Training loop
    print("Starting training...\n")
    best_val_metric = -1
    
    for epoch in range(num_epochs):
        print(f"Epoch {epoch+1}/{num_epochs}")
        print("-" * 80)
        
        # Train
        train_loss = train_epoch(model, train_loader, optimizer, device)
        
        # Validate
        pos_sim, neg_sim = validate(model, val_loader, device)
        
        print(f"\nResults:")
        print(f"  Train loss: {train_loss:.4f}")
        print(f"  Val positive similarity: {pos_sim:.4f}")
        print(f"  Val negative similarity: {neg_sim:.4f}")
        print(f"  Separation: {pos_sim - neg_sim:.4f}")
        
        # Save checkpoint
        separation = pos_sim - neg_sim
        if separation > best_val_metric:
            best_val_metric = separation
            checkpoint_path = checkpoint_dir / 'best_model_10M.pt'
            model.save_checkpoint(
                checkpoint_path,
                metadata={
                    'epoch': epoch + 1,
                    'train_loss': train_loss,
                    'pos_sim': pos_sim,
                    'neg_sim': neg_sim,
                    'separation': separation,
                    'dataset': '10M'
                }
            )
            print(f"  ✅ Saved best model (separation: {best_val_metric:.4f})")
        
        print()
    
    print("="*80)
    print("✅ Training complete!")
    print(f"   Best separation: {best_val_metric:.4f}")
    print("="*80)


if __name__ == "__main__":
    main()
