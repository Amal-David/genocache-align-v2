# Use the same NALEncoder from training
import torch
import torch.nn as nn
import torch.nn.functional as F

class NALEncoder(nn.Module):
    def __init__(self, emb_dim=128, seed_len=512, vocab_size=5, 
                 hidden_dim=128, num_layers=4):
        super().__init__()
        self.emb_dim = emb_dim
        self.seed_len = seed_len
        self.vocab_size = vocab_size
        
        self.embedding = nn.Embedding(vocab_size, hidden_dim)
        self.conv_blocks = nn.ModuleList([
            nn.Conv1d(hidden_dim, hidden_dim, kernel_size=7, padding=3)
            for _ in range(num_layers)
        ])
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Linear(hidden_dim, emb_dim)
        self.norm = nn.LayerNorm(emb_dim)
    
    def forward(self, x):
        x = self.embedding(x)
        x = x.transpose(1, 2)
        
        for conv in self.conv_blocks:
            x = F.relu(conv(x))
        
        x = self.pool(x).squeeze(-1)
        x = self.fc(x)
        x = self.norm(x)
        x = F.normalize(x, p=2, dim=-1)
        return x
