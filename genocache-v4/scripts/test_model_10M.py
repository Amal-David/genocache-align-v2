#!/usr/bin/env python3
"""
Comprehensive validation testing for 10M trained model
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
from pathlib import Path
import h5py
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder

def seq_to_tensor(seq):
    """Convert DNA sequence to tensor"""
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    tensor = torch.zeros(512, dtype=torch.long)
    for i, base in enumerate(seq[:512]):
        tensor[i] = base_to_idx.get(base.upper(), 4)
    return tensor

def test_basic_patterns(model, device):
    """Test basic sequence patterns"""
    print("="*80)
    print("Basic Pattern Tests")
    print("="*80)
    
    tests = [
        ("Identical sequences", "A"*512, "A"*512, ">0.95"),
        ("99% similar", "ACGT"*128, ("ACGT"*127 + "T"*5 + "ACGT"*1), ">0.80"),
        ("Very different", "A"*512, "T"*512, "<0.30"),
    ]
    
    for name, seq1, seq2, expected in tests:
        t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
        t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
        with torch.no_grad():
            e1 = model(t1)
            e2 = model(t2)
            sim = F.cosine_similarity(e1, e2).item()
        status = "✅" if eval(f"{sim} {expected}") else "⚠️"
        print(f"{status} {name:25s}: {sim:7.4f}  (expected {expected})")
    print()

def test_random_sequences(model, device, num_tests=100):
    """Test random unrelated sequences (key test for hard negatives)"""
    print("="*80)
    print(f"Random Sequence Test ({num_tests} pairs)")
    print("="*80)
    
    np.random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    similarities = []
    
    for _ in tqdm(range(num_tests), desc="Testing random pairs"):
        seq1 = ''.join(np.random.choice(bases, 512))
        seq2 = ''.join(np.random.choice(bases, 512))
        
        t1 = seq_to_tensor(seq1).unsqueeze(0).to(device)
        t2 = seq_to_tensor(seq2).unsqueeze(0).to(device)
        
        with torch.no_grad():
            e1 = model(t1)
            e2 = model(t2)
            sim = F.cosine_similarity(e1, e2).item()
        similarities.append(sim)
    
    mean_sim = np.mean(similarities)
    std_sim = np.std(similarities)
    max_sim = np.max(similarities)
    
    status = "✅" if mean_sim < 0.5 else "⚠️"
    print(f"{status} Mean similarity: {mean_sim:.4f} ± {std_sim:.4f}")
    print(f"   Max similarity:  {max_sim:.4f}")
    print(f"   Expected: Mean < 0.5 (lower is better)")
    print()
    
    return mean_sim

def test_validation_set(model, device, data_path, num_samples=1000):
    """Test on validation set from training data"""
    print("="*80)
    print(f"Validation Set Test ({num_samples} examples)")
    print("="*80)
    
    with h5py.File(data_path, 'r') as f:
        total = len(f['anchors'])
        # Use last 1000 as validation
        indices = np.random.choice(range(total-10000, total), num_samples, replace=False)
        
        pos_sims = []
        neg_sims = []
        
        for idx in tqdm(indices, desc="Testing validation set"):
            anchor = f['anchors'][idx].decode('utf-8')
            positive = f['positives'][idx].decode('utf-8')
            negative = f['negatives'][idx].decode('utf-8')
            
            t_anchor = seq_to_tensor(anchor).unsqueeze(0).to(device)
            t_pos = seq_to_tensor(positive).unsqueeze(0).to(device)
            t_neg = seq_to_tensor(negative).unsqueeze(0).to(device)
            
            with torch.no_grad():
                e_anchor = model(t_anchor)
                e_pos = model(t_pos)
                e_neg = model(t_neg)
                
                pos_sim = F.cosine_similarity(e_anchor, e_pos).item()
                neg_sim = F.cosine_similarity(e_anchor, e_neg).item()
                
                pos_sims.append(pos_sim)
                neg_sims.append(neg_sim)
    
    mean_pos = np.mean(pos_sims)
    mean_neg = np.mean(neg_sims)
    separation = mean_pos - mean_neg
    
    print(f"Positive similarity: {mean_pos:.4f} (target: >0.8)")
    print(f"Negative similarity: {mean_neg:.4f} (target: <0.4)")
    print(f"Separation:          {separation:.4f} (target: >0.5)")
    
    if separation > 0.5:
        print("✅ EXCELLENT: Model learned strong separation!")
    elif separation > 0.3:
        print("⚠️ OK: Model learned some separation, could be better")
    else:
        print("❌ POOR: Model did not learn good separation")
    print()
    
    return mean_pos, mean_neg, separation

def test_error_tolerance(model, device):
    """Test tolerance to different error rates"""
    print("="*80)
    print("Error Tolerance Test")
    print("="*80)
    
    np.random.seed(42)
    base_seq = ''.join(np.random.choice(['A','C','G','T'], 512))
    
    error_rates = [0.01, 0.03, 0.05, 0.08, 0.10]
    
    for error_rate in error_rates:
        # Introduce errors
        mutated = list(base_seq)
        num_errors = int(512 * error_rate)
        error_positions = np.random.choice(512, num_errors, replace=False)
        
        for pos in error_positions:
            bases = ['A', 'C', 'G', 'T']
            bases.remove(mutated[pos])
            mutated[pos] = np.random.choice(bases)
        
        mutated_seq = ''.join(mutated)
        
        t1 = seq_to_tensor(base_seq).unsqueeze(0).to(device)
        t2 = seq_to_tensor(mutated_seq).unsqueeze(0).to(device)
        
        with torch.no_grad():
            e1 = model(t1)
            e2 = model(t2)
            sim = F.cosine_similarity(e1, e2).item()
        
        status = "✅" if sim > 0.7 else "⚠️"
        print(f"{status} {error_rate*100:4.1f}% errors: similarity = {sim:.4f}")
    print()

def main():
    print("="*80)
    print("GenoCache V4 - Comprehensive Model Validation")
    print("="*80)
    print()
    
    # Load model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GenoCacheEncoder(emb_dim=256, seed_len=512)
    
    checkpoint_path = 'models/checkpoints/best_model_10M.pt'
    if not Path(checkpoint_path).exists():
        print(f"❌ Checkpoint not found: {checkpoint_path}")
        return
    
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(device)
    model.eval()
    
    print(f"✅ Loaded model from {checkpoint_path}")
    print(f"   Training epoch: {checkpoint['metadata']['epoch']}")
    print(f"   Training separation: {checkpoint['metadata']['separation']:.4f}")
    print()
    
    # Run all tests
    test_basic_patterns(model, device)
    
    random_sim = test_random_sequences(model, device, num_tests=100)
    
    data_path = 'data/training_10M.h5'
    if Path(data_path).exists():
        pos_sim, neg_sim, sep = test_validation_set(model, device, data_path, num_samples=1000)
    
    test_error_tolerance(model, device)
    
    # Final summary
    print("="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    if Path(data_path).exists():
        print(f"Validation set separation: {sep:.4f}")
    print(f"Random sequence similarity: {random_sim:.4f}")
    print()
    
    # Decision criteria
    if Path(data_path).exists() and sep > 0.6 and random_sim < 0.4:
        print("✅ EXCELLENT: Model ready for curriculum learning (130M)")
        print("   - Strong separation on validation set")
        print("   - Good handling of random negatives")
        print("   - Next: Scale to 130M with curriculum (0%→10% errors)")
    elif Path(data_path).exists() and sep > 0.4:
        print("✅ GOOD: Model learned meaningful patterns")
        print("   - Decent separation")
        print("   - May benefit from more training or curriculum")
    else:
        print("⚠️ NEEDS IMPROVEMENT: Consider adjustments")
        print("   - May need different architecture")
        print("   - Or different training strategy")
    
    print()
    print("="*80)

if __name__ == "__main__":
    main()
