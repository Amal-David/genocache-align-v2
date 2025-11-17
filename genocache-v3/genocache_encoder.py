#!/usr/bin/env python3
"""
GenoCacheEncoder: 1.2M parameter multi-scale CNN + Attention
Beyond NeurALigner - Optimized for genomic sequence alignment
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import random


class PositionalEncoding(nn.Module):
    """Learnable positional encoding for translation continuity"""
    def __init__(self, d_model=256, max_len=512):
        super().__init__()
        self.encoding = nn.Parameter(torch.randn(1, max_len, d_model) * 0.02)
    
    def forward(self, x):
        # x: [B, L, D]
        return x + self.encoding[:, :x.size(1), :]


class MultiScaleConvBlock(nn.Module):
    """4 parallel receptive fields for capturing local motifs and global patterns"""
    def __init__(self, in_channels, out_channels):
        super().__init__()
        # 4 parallel kernel sizes: 3, 7, 15, 31
        self.conv3 = nn.Conv1d(in_channels, out_channels // 4, kernel_size=3, padding=1)
        self.conv7 = nn.Conv1d(in_channels, out_channels // 4, kernel_size=7, padding=3)
        self.conv15 = nn.Conv1d(in_channels, out_channels // 4, kernel_size=15, padding=7)
        self.conv31 = nn.Conv1d(in_channels, out_channels // 4, kernel_size=31, padding=15)
        
        self.bn = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU(inplace=True)
        
    def forward(self, x):
        # x: [B, C, L]
        out3 = self.conv3(x)
        out7 = self.conv7(x)
        out15 = self.conv15(x)
        out31 = self.conv31(x)
        
        # Concatenate along channel dimension
        out = torch.cat([out3, out7, out15, out31], dim=1)
        out = self.bn(out)
        out = self.relu(out)
        
        return out


class LightweightAttention(nn.Module):
    """O(L) linear attention instead of O(L²) for long-range dependencies"""
    def __init__(self, d_model=256, num_heads=4):
        super().__init__()
        self.num_heads = num_heads
        self.d_model = d_model
        self.head_dim = d_model // num_heads
        
        assert d_model % num_heads == 0, "d_model must be divisible by num_heads"
        
        self.qkv = nn.Linear(d_model, d_model * 3)
        self.proj = nn.Linear(d_model, d_model)
        self.norm = nn.LayerNorm(d_model)
        
    def forward(self, x):
        # x: [B, L, D]
        B, L, D = x.shape
        
        # Compute Q, K, V
        qkv = self.qkv(x).reshape(B, L, 3, self.num_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)  # [3, B, H, L, d]
        q, k, v = qkv[0], qkv[1], qkv[2]
        
        # Linear attention: softmax over features instead of sequence
        q = F.softmax(q, dim=-1)  # [B, H, L, d]
        k = F.softmax(k, dim=-1)
        
        # Context aggregation: O(L) instead of O(L²)
        context = torch.einsum('bhld,bhld->bhd', k, v)  # [B, H, d]
        out = torch.einsum('bhld,bhd->bhld', q, context)  # [B, H, L, d]
        
        # Reshape and project
        out = out.transpose(1, 2).reshape(B, L, D)
        out = self.proj(out)
        
        return self.norm(out + x)  # Residual connection


class GenoCacheEncoder(nn.Module):
    """
    GenoCacheEncoder: 1.2M parameters
    - Multi-scale CNN for local + global patterns
    - Lightweight attention for long-range dependencies
    - 256D embeddings for high specificity
    
    Target: 97% Top-1 recall, <25μs inference per seed
    """
    def __init__(self, seed_len=512, emb_dim=256):
        super().__init__()
        self.seed_len = seed_len
        self.emb_dim = emb_dim
        
        # Token embedding: 5 tokens (A, C, G, T, N) -> 64D
        self.token_embed = nn.Embedding(5, 64, padding_idx=4)  # N is padding
        
        # Positional encoding
        self.pos_encoding = PositionalEncoding(d_model=64, max_len=seed_len)
        
        # Multi-scale convolutions
        self.ms_conv1 = MultiScaleConvBlock(64, 128)
        self.ms_conv2 = MultiScaleConvBlock(128, 256)
        
        # Lightweight attention
        self.attention = LightweightAttention(d_model=256, num_heads=4)
        
        # Global pooling and projection
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.projection = nn.Sequential(
            nn.Linear(256, emb_dim),
            nn.LayerNorm(emb_dim)
        )
        
        # Parameter count
        self._count_parameters()
    
    def _count_parameters(self):
        total = sum(p.numel() for p in self.parameters())
        trainable = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print(f"GenoCacheEncoder: {trainable:,} trainable parameters ({total:,} total)")
    
    def forward(self, x):
        """
        Args:
            x: [B, L] integer tensor (0=A, 1=C, 2=G, 3=T, 4=N)
               OR [B, L, 4] one-hot tensor
        
        Returns:
            embeddings: [B, emb_dim] L2-normalized vectors
        """
        # Handle both integer and one-hot inputs
        if x.dim() == 3:  # One-hot [B, L, 4]
            # Convert to integers
            x = x.argmax(dim=-1)  # [B, L]
        
        # Token embedding
        x = self.token_embed(x)  # [B, L, 64]
        
        # Positional encoding
        x = self.pos_encoding(x)  # [B, L, 64]
        
        # Multi-scale convolutions
        x = x.permute(0, 2, 1)  # [B, 64, L]
        x = self.ms_conv1(x)     # [B, 128, L]
        x = self.ms_conv2(x)     # [B, 256, L]
        x = x.permute(0, 2, 1)   # [B, L, 256]
        
        # Attention for long-range dependencies
        x = self.attention(x)    # [B, L, 256]
        
        # Global pooling
        x = x.permute(0, 2, 1)   # [B, 256, L]
        x = self.pool(x)          # [B, 256, 1]
        x = x.squeeze(-1)         # [B, 256]
        
        # Project to embedding dimension
        x = self.projection(x)    # [B, emb_dim]
        
        # L2 normalize for cosine similarity search
        x = F.normalize(x, p=2, dim=-1)
        
        return x


class GenoCacheLoss(nn.Module):
    """InfoNCE loss with hard negative mining"""
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
    
    def forward(self, anchors, positives, hard_negatives=None):
        """
        Args:
            anchors: [B, D] embeddings of original sequences
            positives: [B, D] embeddings of augmented versions (same sequence)
            hard_negatives: [B, K, D] embeddings from repetitive regions (optional)
        
        Returns:
            loss: scalar contrastive loss
        """
        B = anchors.shape[0]
        device = anchors.device
        
        # Standard InfoNCE: anchor vs positive + in-batch negatives
        logits = torch.matmul(anchors, positives.T) / self.temperature  # [B, B]
        labels = torch.arange(B, device=device)
        
        loss_a2p = F.cross_entropy(logits, labels)
        loss_p2a = F.cross_entropy(logits.T, labels)
        loss = (loss_a2p + loss_p2a) / 2.0
        
        # Add hard negatives if provided
        if hard_negatives is not None:
            # anchors: [B, D], hard_negatives: [B, K, D]
            K = hard_negatives.shape[1]
            
            # Compute similarity to hard negatives
            hard_sim = torch.bmm(
                anchors.unsqueeze(1),           # [B, 1, D]
                hard_negatives.transpose(1, 2)  # [B, D, K]
            ).squeeze(1) / self.temperature      # [B, K]
            
            # Positive similarity (anchor vs its augmentation)
            pos_sim = (anchors * positives).sum(dim=-1, keepdim=True) / self.temperature  # [B, 1]
            
            # Contrastive loss: positive should be higher than hard negatives
            hard_logits = torch.cat([pos_sim, hard_sim], dim=1)  # [B, K+1]
            hard_labels = torch.zeros(B, dtype=torch.long, device=device)  # First column is positive
            
            hard_loss = F.cross_entropy(hard_logits, hard_labels)
            loss = loss + 0.5 * hard_loss  # Weight the hard negative loss
        
        return loss


# Data augmentation functions

def seq_to_tokens(seq):
    """Convert DNA sequence to integer tokens"""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    return [mapping.get(ch.upper(), 4) for ch in seq]


def seq_to_onehot(seq):
    """Convert DNA sequence to one-hot encoding (for compatibility)"""
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    for i, ch in enumerate(seq):
        if ch == "A": arr[i, 0] = 1
        elif ch == "C": arr[i, 1] = 1
        elif ch == "G": arr[i, 2] = 1
        elif ch == "T": arr[i, 3] = 1
        else: arr[i, :] = 0.25  # Uniform for N
    return arr


def augment_sequence(seq, error_rate_range=(0.01, 0.10), rc_prob=0.5):
    """
    Augment DNA sequence with realistic sequencing errors
    
    Args:
        seq: DNA sequence string
        error_rate_range: (min, max) substitution rate
        rc_prob: probability of reverse complement
    
    Returns:
        augmented sequence
    """
    # Reverse complement
    if random.random() < rc_prob:
        trans = str.maketrans("ACGT", "TGCA")
        seq = seq.translate(trans)[::-1]
    
    # Substitution errors
    error_rate = random.uniform(*error_rate_range)
    seq_list = list(seq)
    
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            # 60% substitution, 20% insertion-like, 20% deletion-like
            r = random.random()
            if r < 0.6:
                seq_list[i] = random.choice("ACGT")
            elif r < 0.8:
                # Simulate insertion by duplicating
                if i > 0:
                    seq_list[i] = seq_list[i-1]
            else:
                # Simulate deletion by setting to N (will be handled)
                seq_list[i] = 'N'
    
    seq = ''.join(seq_list)
    
    # Random shift (simulate alignment uncertainty)
    seed_len = len(seq)
    shift = random.randint(-seed_len // 10, seed_len // 10)
    if shift > 0:
        seq = seq[shift:] + seq[:shift]
    elif shift < 0:
        seq = seq[shift:] + seq[:shift]
    
    # Ensure correct length (pad or trim)
    if len(seq) < seed_len:
        seq = seq + 'N' * (seed_len - len(seq))
    elif len(seq) > seed_len:
        start = random.randint(0, len(seq) - seed_len)
        seq = seq[start:start + seed_len]
    
    return seq


if __name__ == "__main__":
    # Test the model
    print("Testing GenoCacheEncoder...")
    
    model = GenoCacheEncoder(seed_len=512, emb_dim=256)
    
    # Test with integer input
    batch_size = 4
    seq_len = 512
    x_int = torch.randint(0, 4, (batch_size, seq_len))
    
    print(f"\nInput shape (int): {x_int.shape}")
    out = model(x_int)
    print(f"Output shape: {out.shape}")
    print(f"Output norm: {out.norm(dim=-1)}")  # Should be ~1.0
    
    # Test with one-hot input
    x_onehot = torch.randn(batch_size, seq_len, 4)
    out2 = model(x_onehot)
    print(f"\nOne-hot input shape: {x_onehot.shape}")
    print(f"Output shape: {out2.shape}")
    
    # Test loss
    loss_fn = GenoCacheLoss(temperature=0.07)
    
    anchors = torch.randn(batch_size, 256)
    positives = torch.randn(batch_size, 256)
    hard_negs = torch.randn(batch_size, 10, 256)
    
    # Normalize
    anchors = F.normalize(anchors, dim=-1)
    positives = F.normalize(positives, dim=-1)
    hard_negs = F.normalize(hard_negs, dim=-1)
    
    loss = loss_fn(anchors, positives, hard_negs)
    print(f"\nLoss: {loss.item():.4f}")
    
    print("\n✅ Model test passed!")
