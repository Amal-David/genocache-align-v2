#!/usr/bin/env python3
"""
NAL-Aligned Encoder - Exact Match to NeuralAligner Design

Key differences from our old encoder:
- 128D output (not 256D!)
- NO positional encoding (NAL uses shift augmentation instead)
- NO attention layers (NAL uses simple conv)
- NO multi-scale complexity
- Simple bidirectional convolutions only

Reference: NeuralAligner paper Section 3.1, A.2
Architecture: Hyena-DNA tiny-1k inspired, 0.5M params
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class NALEncoder(nn.Module):
    """
    NAL-aligned encoder - matches NeuralAligner's design exactly
    
    Architecture:
    - Token embedding (vocab_size=5 → hidden_dim)
    - 4 bidirectional conv layers (kernel=7)
    - Average pooling over sequence length
    - L2 normalization
    - Optional projection head (training only)
    
    Output: 128D embeddings (matches NAL)
    Parameters: ~0.5M (matches NAL)
    """
    
    def __init__(
        self,
        emb_dim=128,        # Output dimension (matches NAL!)
        seed_len=512,       # Max sequence length
        vocab_size=5,       # A, C, G, T, N
        hidden_dim=128,     # Hidden dimension (single scale)
        num_layers=4,       # Depth of conv stack
        kernel_size=7,      # Conv kernel (bidirectional)
        dropout=0.1
    ):
        super().__init__()
        
        self.emb_dim = emb_dim
        self.seed_len = seed_len
        self.hidden_dim = hidden_dim
        
        # Token embedding (vocab → hidden)
        self.token_embedding = nn.Embedding(vocab_size, hidden_dim)
        
        # Bidirectional conv layers (NAL design)
        # Padding = kernel_size//2 makes it bidirectional
        self.conv_layers = nn.ModuleList()
        for i in range(num_layers):
            self.conv_layers.append(nn.Sequential(
                nn.Conv1d(
                    hidden_dim, 
                    hidden_dim, 
                    kernel_size=kernel_size, 
                    padding=kernel_size//2,  # Bidirectional!
                    bias=False
                ),
                nn.BatchNorm1d(hidden_dim),
                nn.GELU(),
                nn.Dropout(dropout)
            ))
        
        # Projection head for contrastive learning (TRAINING ONLY)
        # Used during training to prevent dimensional collapse
        # NOT used during inference/indexing!
        self.projection = nn.Sequential(
            nn.Linear(hidden_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )
        
        self._init_weights()
    
    def _init_weights(self):
        """Initialize weights properly"""
        for m in self.modules():
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
    
    def forward(self, x, use_projection=False):
        """
        Forward pass
        
        Args:
            x: (B, L) tensor of token indices
            use_projection: If True, apply projection head (training only!)
        
        Returns:
            embeddings: (B, emb_dim) tensor
        """
        B, L = x.shape
        
        # Token embedding: (B, L) → (B, L, hidden_dim)
        x = self.token_embedding(x)
        
        # Transpose for conv: (B, L, hidden_dim) → (B, hidden_dim, L)
        x = x.transpose(1, 2)
        
        # Bidirectional conv layers
        for conv_layer in self.conv_layers:
            x = conv_layer(x) + x  # Residual connection
        
        # Average pooling over sequence: (B, hidden_dim, L) → (B, hidden_dim)
        x = torch.mean(x, dim=2)
        
        # L2 normalization (NAL does this!)
        x = F.normalize(x, p=2, dim=1)
        
        # Apply projection head if requested (training only!)
        if use_projection:
            x = self.projection(x)
            # Normalize again after projection
            x = F.normalize(x, p=2, dim=1)
        
        return x
    
    def encode_sequences(self, sequences, device='cuda', use_projection=False):
        """
        Encode multiple DNA sequences
        
        Args:
            sequences: List of DNA strings
            device: Device to use
            use_projection: Whether to use projection head
        
        Returns:
            embeddings: (N, emb_dim) tensor
        """
        self.eval()
        
        # Convert sequences to tensors
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        
        tensors = []
        for seq in sequences:
            indices = [base_to_idx.get(b, 4) for b in seq]
            # Pad or trim to seed_len
            if len(indices) < self.seed_len:
                indices = indices + [4] * (self.seed_len - len(indices))
            else:
                indices = indices[:self.seed_len]
            tensors.append(torch.tensor(indices, dtype=torch.long))
        
        x = torch.stack(tensors).to(device)
        
        with torch.no_grad():
            embeddings = self.forward(x, use_projection=use_projection)
        
        return embeddings


class InfoNCELoss(nn.Module):
    """
    InfoNCE loss for contrastive learning (exact NAL formula)
    
    Reference: NeuralAligner Equation 1
    """
    
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
    
    def forward(self, anchor_emb, positive_emb):
        """
        Compute InfoNCE loss
        
        Args:
            anchor_emb: (B, D) embeddings (already normalized)
            positive_emb: (B, D) embeddings (already normalized)
        
        Returns:
            loss: scalar
            pos_sim: average positive similarity
            neg_sim: average negative similarity
        """
        batch_size = anchor_emb.shape[0]
        
        # Compute similarity matrix: (B, B)
        # similarity[i,j] = anchor[i] · positive[j] / temp
        similarity = torch.mm(anchor_emb, positive_emb.T) / self.temperature
        
        # Positive pairs: diagonal elements
        pos_sim = torch.diagonal(similarity)
        
        # InfoNCE loss: -log(exp(pos) / sum_exp(all))
        log_sum_exp = torch.logsumexp(similarity, dim=1)
        loss = -torch.mean(pos_sim - log_sum_exp)
        
        # Compute average similarities for monitoring
        pos_sim_avg = torch.mean(pos_sim).item()
        
        # Negative similarities (off-diagonal)
        mask = torch.eye(batch_size, device=similarity.device).bool()
        neg_similarities = similarity.masked_fill(mask, float('-inf'))
        neg_sim_avg = torch.mean(neg_similarities[~mask]).item()
        
        return loss, pos_sim_avg, neg_sim_avg


def test_encoder():
    """Test NAL encoder"""
    print("Testing NAL-Aligned Encoder")
    print("=" * 60)
    
    # Create model
    model = NALEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dim=128,
        num_layers=4
    )
    
    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
    print(f"Output dimension: {model.emb_dim}")
    
    # Test forward pass
    batch_size = 16
    seq_len = 512
    x = torch.randint(0, 5, (batch_size, seq_len))
    
    # Without projection (inference/indexing)
    emb_no_proj = model(x, use_projection=False)
    print(f"\nWithout projection: {emb_no_proj.shape}")
    print(f"  Mean norm: {torch.norm(emb_no_proj, dim=1).mean():.4f} (should be ~1.0)")
    
    # With projection (training)
    emb_with_proj = model(x, use_projection=True)
    print(f"With projection: {emb_with_proj.shape}")
    print(f"  Mean norm: {torch.norm(emb_with_proj, dim=1).mean():.4f} (should be ~1.0)")
    
    # Test InfoNCE loss
    criterion = InfoNCELoss(temperature=0.07)
    loss, pos_sim, neg_sim = criterion(emb_with_proj, emb_with_proj)
    print(f"\nInfoNCE loss test:")
    print(f"  Loss: {loss:.4f}")
    print(f"  Pos sim: {pos_sim:.4f}")
    print(f"  Neg sim: {neg_sim:.4f}")
    print(f"  Separation: {pos_sim - neg_sim:.4f}")
    
    print("\n✅ NAL encoder test complete!")


if __name__ == '__main__':
    test_encoder()
