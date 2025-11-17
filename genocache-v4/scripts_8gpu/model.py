"""
GenoCache V4 - Model Architecture
Hyena-DNA inspired encoder for DNA sequence embeddings
"""

import torch
import torch.nn as nn
import math

class HyenaConv1d(nn.Module):
    """1D Convolution with GELU activation"""
    def __init__(self, in_channels, out_channels, kernel_size=3):
        super().__init__()
        self.conv = nn.Conv1d(
            in_channels, 
            out_channels, 
            kernel_size,
            padding=kernel_size // 2
        )
        self.activation = nn.GELU()
    
    def forward(self, x):
        return self.activation(self.conv(x))

class HyenaBlock(nn.Module):
    """Hyena block with convolution and feedforward"""
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
    """
    GenoCache V4 Model
    
    Architecture:
    - Embedding layer (vocab=5 for A,C,G,T,N)
    - Positional embedding
    - 6 Hyena blocks
    - Global average pooling
    - L2 normalization
    
    Output: 128D L2-normalized embeddings
    """
    def __init__(self, vocab_size=5, embed_dim=128, num_layers=6, max_len=512):
        super().__init__()
        
        # Token embedding
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=4)
        
        # Positional embedding (learnable)
        self.pos_embedding = nn.Parameter(torch.randn(1, max_len, embed_dim) * 0.02)
        
        # Hyena blocks
        self.blocks = nn.ModuleList([
            HyenaBlock(embed_dim, kernel_size=3)
            for _ in range(num_layers)
        ])
        
        # Final norm
        self.norm = nn.LayerNorm(embed_dim)
        
        # Initialize weights
        self.apply(self._init_weights)
    
    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.02)
            if module.bias is not None:
                torch.nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.02)
        elif isinstance(module, nn.LayerNorm):
            torch.nn.init.zeros_(module.bias)
            torch.nn.init.ones_(module.weight)
    
    def forward(self, x):
        """
        Args:
            x: (batch, seq_len) - Token indices
        
        Returns:
            embeddings: (batch, embed_dim) - L2-normalized embeddings
        """
        batch_size, seq_len = x.shape
        
        # Token embedding
        x = self.embedding(x)  # (batch, seq_len, embed_dim)
        
        # Add positional embedding
        x = x + self.pos_embedding[:, :seq_len, :]
        
        # Process through Hyena blocks
        for block in self.blocks:
            x = block(x)
        
        # Final norm
        x = self.norm(x)
        
        # Global average pooling
        x = x.mean(dim=1)  # (batch, embed_dim)
        
        # L2 normalization
        x = torch.nn.functional.normalize(x, p=2, dim=1)
        
        return x
    
    def get_num_params(self):
        """Count total parameters"""
        return sum(p.numel() for p in self.parameters())

def load_model(checkpoint_path, device='cuda'):
    """Load model from checkpoint"""
    model = GenoCache(
        vocab_size=5,
        embed_dim=128,
        num_layers=6,
        max_len=512
    )
    
    checkpoint = torch.load(checkpoint_path, map_location=device)
    
    # Handle different checkpoint formats
    if 'model_state_dict' in checkpoint:
        state_dict = checkpoint['model_state_dict']
    else:
        state_dict = checkpoint
    
    # Remove 'module.' prefix if present (from DDP)
    new_state_dict = {}
    for k, v in state_dict.items():
        if k.startswith('module.'):
            new_state_dict[k[7:]] = v
        else:
            new_state_dict[k] = v
    
    model.load_state_dict(new_state_dict)
    model = model.to(device)
    model.eval()
    
    return model

if __name__ == '__main__':
    # Test model
    model = GenoCache()
    print(f"Model parameters: {model.get_num_params():,}")
    
    # Test forward pass
    x = torch.randint(0, 5, (4, 512))
    with torch.no_grad():
        embeddings = model(x)
    print(f"Input shape: {x.shape}")
    print(f"Output shape: {embeddings.shape}")
    print(f"Output norm: {embeddings.norm(dim=1)}")  # Should be ~1.0
