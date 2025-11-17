#!/usr/bin/env python3
"""
Test the trained model to verify embeddings are meaningful
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder

def seq_to_tensor(seq):
    """Convert DNA sequence to tensor"""
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    tensor = torch.zeros(512, dtype=torch.long)
    for i, base in enumerate(seq[:512]):
        tensor[i] = base_to_idx.get(base.upper(), 4)
    return tensor

def test_model():
    print("="*80)
    print("Testing Trained Model")
    print("="*80)
    
    # Load model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GenoCacheEncoder(emb_dim=256, seed_len=512)
    
    checkpoint_path = 'models/checkpoints/best_model.pt'
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(device)
    model.eval()
    
    print(f"✅ Loaded model from {checkpoint_path}")
    print(f"   Epoch: {checkpoint['metadata']['epoch']}")
    print(f"   Separation: {checkpoint['metadata']['separation']:.4f}")
    print()
    
    # Test sequences
    print("Test 1: Identical sequences")
    seq1 = "A" * 512
    seq2 = "A" * 512
    t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
    t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
    with torch.no_grad():
        e1 = model(t1)
        e2 = model(t2)
        sim = F.cosine_similarity(e1, e2).item()
    print(f"  Similarity: {sim:.4f} (should be ~1.0)")
    
    # Test 2: Slightly different (1% errors)
    print("\nTest 2: Similar sequences (99% identical)")
    seq1 = "ACGTACGT" * 64
    seq2 = list(seq1)
    for i in range(0, 512, 100):
        seq2[i] = 'T'  # Change ~1%
    seq2 = ''.join(seq2)
    t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
    t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
    with torch.no_grad():
        e1 = model(t1)
        e2 = model(t2)
        sim = F.cosine_similarity(e1, e2).item()
    print(f"  Similarity: {sim:.4f} (should be >0.8)")
    
    # Test 3: Very different
    print("\nTest 3: Different sequences")
    seq1 = "A" * 512
    seq2 = "T" * 512
    t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
    t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
    with torch.no_grad():
        e1 = model(t1)
        e2 = model(t2)
        sim = F.cosine_similarity(e1, e2).item()
    print(f"  Similarity: {sim:.4f} (should be <0.5)")
    
    # Test 4: Random sequences
    print("\nTest 4: Random unrelated sequences")
    np.random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    seq1 = ''.join(np.random.choice(bases, 512))
    seq2 = ''.join(np.random.choice(bases, 512))
    t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
    t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
    with torch.no_grad():
        e1 = model(t1)
        e2 = model(t2)
        sim = F.cosine_similarity(e1, e2).item()
    print(f"  Similarity: {sim:.4f} (should be <0.5)")
    
    # Test 5: Reverse complement
    print("\nTest 5: Reverse complement")
    seq1 = "ACGTACGT" * 64
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq2 = ''.join([complement[b] for b in seq1[::-1]])
    t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
    t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
    with torch.no_grad():
        e1 = model(t1)
        e2 = model(t2)
        sim = F.cosine_similarity(e1, e2).item()
    print(f"  Similarity: {sim:.4f} (with RC augmentation, could be variable)")
    
    print("\n" + "="*80)
    print("✅ Model testing complete!")
    print("="*80)

if __name__ == "__main__":
    test_model()
