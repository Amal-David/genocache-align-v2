"""
GenoCache-Align V2: Fast CNN Tier
Optimized for unique mappers (~70-80% of reads)

Architecture: ImprovedCNN
- Target: 95%+ recall on unique regions
- Speed: 1M+ reads/sec with GPU batching
- Index: <50MB compressed per genome
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import json
from pathlib import Path
from datetime import datetime

# Configuration
EMB_DIM = 256  # Embedding dimension (balance between accuracy and index size)
SEED_LEN = 512  # Seed length (bp)


class ImprovedCNN(nn.Module):
    """
    Fast CNN encoder for Tier 1 alignment
    
    Architecture improvements over baseline:
    - Deeper network (5 conv layers vs 3)
    - BatchNorm for training stability
    - Residual connections for gradient flow
    - Adaptive pooling for consistent output
    
    Performance targets:
    - Top-1 accuracy: >64% (proven in h100 validation)
    - Top-10 accuracy: >85%
    - Inference: <1ms per read on GPU
    """
    
    def __init__(self, out_dim=EMB_DIM, input_len=SEED_LEN):
        super().__init__()
        
        self.input_len = input_len
        self.out_dim = out_dim
        
        # Stage 1: Initial feature extraction (4 → 64 channels)
        self.conv1 = nn.Conv1d(4, 64, kernel_size=11, padding=5)
        self.bn1 = nn.BatchNorm1d(64)
        
        # Stage 2: Feature expansion (64 → 128 channels)
        self.conv2 = nn.Conv1d(64, 128, kernel_size=9, padding=4)
        self.bn2 = nn.BatchNorm1d(128)
        
        # Stage 3: Deep features (128 → 256 channels)
        self.conv3 = nn.Conv1d(128, 256, kernel_size=7, padding=3)
        self.bn3 = nn.BatchNorm1d(256)
        
        # Stage 4: Feature refinement (256 → 256 channels, with residual)
        self.conv4 = nn.Conv1d(256, 256, kernel_size=5, padding=2)
        self.bn4 = nn.BatchNorm1d(256)
        
        # Stage 5: Final features (256 → 320 channels)
        self.conv5 = nn.Conv1d(256, 320, kernel_size=3, padding=1)
        self.bn5 = nn.BatchNorm1d(320)
        
        # Global pooling to fixed-size representation
        self.adaptive_pool = nn.AdaptiveAvgPool1d(1)
        
        # Projection to embedding space
        self.fc = nn.Linear(320, out_dim)
        
        # Dropout for regularization
        self.dropout = nn.Dropout(0.1)
        
    def forward(self, x):
        """
        Args:
            x: (batch, 4, seq_len) - One-hot encoded DNA sequences
        Returns:
            embeddings: (batch, out_dim) - L2-normalized embeddings
        """
        # Stage 1
        x = F.relu(self.bn1(self.conv1(x)))
        x = self.dropout(x)
        
        # Stage 2
        x = F.relu(self.bn2(self.conv2(x)))
        x = self.dropout(x)
        
        # Stage 3
        x = F.relu(self.bn3(self.conv3(x)))
        x = self.dropout(x)
        
        # Stage 4 (with residual connection)
        identity = x
        x = F.relu(self.bn4(self.conv4(x)))
        x = x + identity  # Residual connection
        x = self.dropout(x)
        
        # Stage 5
        x = F.relu(self.bn5(self.conv5(x)))
        
        # Global pooling
        x = self.adaptive_pool(x)  # (batch, 320, 1)
        x = x.squeeze(-1)  # (batch, 320)
        
        # Project to embedding space
        x = self.fc(x)  # (batch, out_dim)
        
        # L2 normalize for cosine similarity
        x = F.normalize(x, p=2, dim=1)
        
        return x
    
    def encode_batch(self, sequences):
        """
        Convenience method for encoding a batch of sequences
        
        Args:
            sequences: List of DNA strings or tensor (batch, 4, seq_len)
        Returns:
            embeddings: (batch, out_dim) numpy array
        """
        if isinstance(sequences, list):
            sequences = self.sequences_to_tensor(sequences)
        
        self.eval()
        with torch.no_grad():
            embeddings = self.forward(sequences)
        
        return embeddings.cpu().numpy()
    
    @staticmethod
    def sequences_to_tensor(sequences):
        """Convert DNA sequences to one-hot encoded tensor"""
        from genomic_utils import one_hot_encode
        tensors = [one_hot_encode(seq) for seq in sequences]
        return torch.stack(tensors)
    
    def save_checkpoint(self, filepath, metadata=None):
        """Save model with metadata for provenance tracking"""
        checkpoint = {
            'model_state_dict': self.state_dict(),
            'architecture': 'ImprovedCNN',
            'emb_dim': self.out_dim,
            'input_len': self.input_len,
            'timestamp': datetime.now().isoformat(),
            'metadata': metadata or {}
        }
        torch.save(checkpoint, filepath)
        
        # Also save manifest
        manifest_path = Path(filepath).with_suffix('.json')
        with open(manifest_path, 'w') as f:
            json.dump({
                'checkpoint_path': str(filepath),
                'architecture': checkpoint['architecture'],
                'emb_dim': checkpoint['emb_dim'],
                'input_len': checkpoint['input_len'],
                'timestamp': checkpoint['timestamp'],
                'metadata': checkpoint['metadata']
            }, f, indent=2)
        
        print(f"✅ Saved checkpoint: {filepath}")
        print(f"✅ Saved manifest: {manifest_path}")
    
    @classmethod
    def load_checkpoint(cls, filepath, device='cuda'):
        """Load model from checkpoint"""
        checkpoint = torch.load(filepath, map_location=device)
        
        model = cls(
            out_dim=checkpoint['emb_dim'],
            input_len=checkpoint['input_len']
        )
        model.load_state_dict(checkpoint['model_state_dict'])
        model.to(device)
        model.eval()
        
        print(f"✅ Loaded {checkpoint['architecture']} from {filepath}")
        print(f"   Timestamp: {checkpoint['timestamp']}")
        print(f"   Embedding dim: {checkpoint['emb_dim']}")
        
        return model, checkpoint['metadata']


class ContrastiveLoss(nn.Module):
    """
    Contrastive loss for training seed embeddings
    Positive pairs: Seeds from same genomic position
    Negative pairs: Seeds from different positions
    """
    
    def __init__(self, temperature=0.07, margin=0.5):
        super().__init__()
        self.temperature = temperature
        self.margin = margin
    
    def forward(self, embeddings, labels):
        """
        Args:
            embeddings: (batch, emb_dim) - L2-normalized
            labels: (batch,) - Position indices (seeds from same position have same label)
        """
        # Compute pairwise cosine similarity
        similarity_matrix = torch.matmul(embeddings, embeddings.T)
        
        # Mask for positive pairs (same label, excluding self)
        labels = labels.unsqueeze(1)
        mask_positive = (labels == labels.T).float()
        mask_positive.fill_diagonal_(0)
        
        # Mask for negative pairs
        mask_negative = 1 - mask_positive
        mask_negative.fill_diagonal_(0)
        
        # Positive loss: maximize similarity for same positions
        pos_similarity = similarity_matrix * mask_positive
        pos_loss = -torch.log(torch.exp(pos_similarity / self.temperature) + 1e-8)
        pos_loss = (pos_loss * mask_positive).sum() / (mask_positive.sum() + 1e-8)
        
        # Negative loss: margin-based separation
        neg_similarity = similarity_matrix * mask_negative
        neg_loss = torch.clamp(neg_similarity - self.margin, min=0)
        neg_loss = (neg_loss * mask_negative).sum() / (mask_negative.sum() + 1e-8)
        
        return pos_loss + neg_loss


# Model information for external use
MODEL_INFO = {
    'name': 'ImprovedCNN',
    'tier': 1,
    'purpose': 'Fast alignment for unique mappers',
    'embedding_dim': EMB_DIM,
    'seed_length': SEED_LEN,
    'target_throughput': '1M+ reads/sec',
    'target_recall': '95%+ on unique regions',
    'rescue_trigger': {
        'top2_ratio': 1.2,  # If top-2 scores are too close
        'min_score': 0.7,   # If max score below threshold
        'max_hits': 5       # If too many high-scoring hits
    }
}


if __name__ == "__main__":
    # Test model
    print("Testing ImprovedCNN architecture...")
    
    model = ImprovedCNN()
    print(f"\nModel: {MODEL_INFO['name']}")
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")
    print(f"Embedding dim: {MODEL_INFO['embedding_dim']}")
    
    # Test forward pass
    batch_size = 32
    x = torch.randn(batch_size, 4, SEED_LEN)
    
    with torch.no_grad():
        embeddings = model(x)
    
    print(f"\nTest forward pass:")
    print(f"  Input: {x.shape}")
    print(f"  Output: {embeddings.shape}")
    print(f"  L2 norm: {torch.norm(embeddings, dim=1).mean():.4f} (should be ~1.0)")
    
    print("\n✅ Architecture ready!")
    print("\nNext steps:")
    print("1. Train with train_fast_cnn.py")
    print("2. Encode reference with encode_reference.py")
    print("3. Build FAISS index with build_index.py")
