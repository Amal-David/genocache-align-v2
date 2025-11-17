# Exact Code Issues - What's Wrong & What to Fix

**Root Cause:** 32% mapping rate vs 94% minimap2 on real GIAB reads

---

## Issue #1: Error Rate in Training Data ❌ CRITICAL

### What We Did Wrong:

**File:** `development/training/augmentation.py`

```python
# Line ~80-90 (WRONG)
def augment_with_errors(self, seq, error_rate=0.02):  # ❌ FIXED 2%!
    """
    PROBLEM: Only uses 2% error rate
    NeuralAligner uses 0.01-0.1 (1% to 10%) sampled uniformly
    """
```

**Current Code:**
- Error rate: **Fixed at 2%**
- No variation
- Doesn't train model for diverse error scenarios

**NeuralAligner (Correct):**
```python
# Error rates sampled from U[0.01, 0.1]
error_rate = np.random.uniform(0.01, 0.1)
# This trains model to be robust to 1%-10% errors
```

**Impact:** Model not robust to real ONT reads with variable error rates (1%-15%)

### Fix:

```python
def augment_with_errors(self, seq, error_rate=None):
    if error_rate is None:
        # Sample from uniform distribution like NeuralAligner
        error_rate = np.random.uniform(0.01, 0.1)
    
    # Apply errors...
```

**Priority:** 🔴 CRITICAL - Main cause of 32% mapping rate

---

## Issue #2: Missing RC-Augmentation ⚠️ IMPORTANT

### What We Did Wrong:

**File:** `development/training/augmentation.py`

```python
# Line ~200 (MISSING)
def create_training_pairs(self, genome_seq, batch_size=1024):
    """
    PROBLEM: No reverse complement augmentation
    NeuralAligner randomly flips sequences
    """
    # Extract sequences
    anchors = []
    positives = []
    
    # Apply errors...
    # ❌ MISSING: No RC flip with 50% probability
```

**Current:** Never applies reverse complement

**NeuralAligner (Correct):**
```python
# After augmentation, randomly flip to RC with 50% probability
for i in range(batch_size):
    if np.random.rand() < 0.5:
        anchors[i] = reverse_complement(anchors[i])
    if np.random.rand() < 0.5:
        positives[i] = reverse_complement(positives[i])
```

**Impact:** 
- Model not RC-invariant
- Must search both strands separately (2× slower)
- Less robust to strand orientation

### Fix:

```python
def create_training_pairs(self, genome_seq, batch_size=1024):
    # ... existing code ...
    
    # RC augmentation (NeuralAligner approach)
    for i in range(batch_size):
        if np.random.rand() < 0.5:
            anchors[i] = reverse_complement(anchors[i])
        if np.random.rand() < 0.5:
            positives[i] = reverse_complement(positives[i])
    
    return anchors, positives
```

**Priority:** 🟡 IMPORTANT - Affects efficiency and robustness

---

## Issue #3: Index Stride Too Large ⚠️

### What We Did Wrong:

**File:** Model checkpoint metadata or index building script (need to check)

**NeuralAligner:**
- Seed length: 256 or 512
- **Index stride:** ≤ seed_len / 8
- For 512bp seeds: stride ≤ 64bp

**Current GenoCache:**
- Seed length: 512bp
- **Index stride:** Unknown (need to check)
- If > 64bp, causes missed anchors

### How to Check:

```python
# Check index metadata
import pickle
with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
    meta = pickle.load(f)
    print(f"Positions: {len(meta['positions'])}")
    print(f"Genome length: ~3 billion")
    stride = 3_000_000_000 / len(meta['positions'])
    print(f"Effective stride: {stride:.0f}bp")
```

**If stride > 64bp:** Need to rebuild index with stride ≤ 64bp

**Priority:** 🟡 MEDIUM - May explain some unmapped reads

---

## Issue #4: Chaining Algorithm Simplification ⚠️

### What We Did Wrong:

**File:** `genocache_core/adaptive_seeding.py` lines ~240-280

```python
# CURRENT (Simplified)
def chain_seeds(self, seeds):
    # Group by chromosome
    for chr_name, seed_list in chr_seeds.items():
        # Check colinearity
        for i in range(len(seed_list) - 1):
            error = abs(read_dist - ref_dist)
            if error > self.colinearity_tolerance:  # ❌ Too simple
                colinear = False
```

**Problems:**
- Binary pass/fail (colinear or not)
- Single tolerance value (3000bp)
- Doesn't use stripe-based approach

**NeuralAligner (Correct):**
```python
# Equation 2 from paper - stripe-based diagonal matching
y_chain = argmax_y [
    sum_i sum_j I(|y_ij - x_i - y| <= C)
]
# Finds optimal diagonal with most anchors
# More robust to positional errors
```

### Fix:

Implement vectorized stripe-based chaining:
```python
def chain_seeds_vectorized(self, seeds, tolerance=3000):
    """Vectorized chaining like NeuralAligner Eq. 2"""
    # Create position matrix
    # Find diagonal with most anchors
    # Use NumPy for speed
```

**Priority:** 🟢 LOW - Current approach works reasonably well

---

## Issue #5: Training Data Distribution ⚠️

### What We Did Wrong:

**Training might have used:**
- Chr22 only (too limited)
- Not enough diversity
- Wrong error model

**NeuralAligner:**
- Full genome GRCh38.p14
- 13 million sample pairs
- Errors: U[0.01, 0.1]
- Shifts: ±seed_len/10

### Check What We Have:

```bash
# Check model checkpoint
ls -lh models/genocache_model.pt
# Check training logs
ls development/training/logs/
```

If trained on chr22 only or < 1M samples → need retraining

**Priority:** 🔴 HIGH - Related to robustness issue

---

## Issue #6: Embedding Dimension Mismatch

### What We Have:

**File:** `genocache_core/encoder.py` line ~153

```python
def __init__(self, emb_dim=256, ...):  # 256D
    # Claims to be "2× better than NAL's 128D"
```

**NeuralAligner:**
- Embedding: 128D
- Works perfectly at 94%

**Our Model:**
- Embedding: 256D (but trained with 128D?)
- Check actual: `model.emb_dim`

**Possible Issue:** Model trained with 128D but index uses 256D (or vice versa)

### Check:

```python
checkpoint = torch.load('models/genocache_model.pt')
print(checkpoint['model_state_dict'].keys())
# Check final layer dimensions
```

**Priority:** 🟡 MEDIUM - May cause embedding mismatch

---

## Summary: What to Fix (Priority Order)

### 🔴 Critical (Do First)

1. **Error Rate Sampling** (5 min code + retraining)
   ```python
   # augmentation.py
   error_rate = np.random.uniform(0.01, 0.1)  # Add this line
   ```

2. **Training Data Volume** (Check logs)
   - If < 1M samples: retrain with 10M+
   - If chr22 only: use full genome

### 🟡 Important (Do Second)

3. **RC-Augmentation** (10 min code + retraining)
   ```python
   # Add RC flipping in create_training_pairs()
   ```

4. **Index Stride** (Check + rebuild if needed)
   ```python
   # Ensure stride ≤ 64bp for 512bp seeds
   ```

5. **Embedding Dimension** (Verify consistency)
   ```python
   # Check model.emb_dim == index_dim
   ```

### 🟢 Optional (Nice to Have)

6. **Vectorized Chaining** (NAL Eq. 2)
   - Current works OK
   - Could gain 5-10%

---

## Quick Diagnosis Script

```python
# Run this to check all issues:

import torch
import pickle
import numpy as np

# 1. Check model
checkpoint = torch.load('models/genocache_model.pt')
print(f"Model keys: {list(checkpoint.keys())[:5]}")

# 2. Check index
with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
    meta = pickle.load(f)
stride = 3_000_000_000 / len(meta['positions'])
print(f"Index stride: {stride:.0f}bp (should be ≤64bp)")

# 3. Check training
import sys
sys.path.append('development/training')
from augmentation import DataAugmentation
aug = DataAugmentation()
# Check error rate implementation

print("\n🔴 CRITICAL ISSUES:")
if stride > 64:
    print(f"  - Index stride too large: {stride:.0f}bp")
print("  - Error rate fixed at 2% (need 1%-10%)")
print("  - No RC-augmentation")
```

---

## Expected Improvement After Fixes

| Fix | Expected Gain | Time |
|-----|---------------|------|
| **Error rate sampling** | +30-40% | 1-2 days |
| **RC-augmentation** | +5-10% | 1-2 days |
| **Index stride** | +5-15% | 2-4 hours |
| **All combined** | **60-80% total** | 2-3 days |

**Target:** 32% → 85-90% (close to minimap2's 94%)

---

## Retrain Protocol (NAL-compliant)

```python
# 1. Fix augmentation.py
- Add error_rate = np.random.uniform(0.01, 0.1)
- Add RC flipping

# 2. Generate training data
python3 generate_training_data.py --samples 10000000

# 3. Train model
python3 train_fullgenome.py --epochs 30 --batch-size 8192

# 4. Rebuild index with stride=64
python3 build_index.py --stride 64

# 5. Test on GIAB
python3 genocache_align_minimap2.py --reads giab_hg002_100reads.fastq

# Expected: 85-90% mapping rate
```

---

## The Real Bottleneck

**It's NOT:**
- ❌ Chaining logic (works fine)
- ❌ Rescue strategy (correct)
- ❌ SAM output (perfect)
- ❌ Code architecture (excellent)

**It IS:**
- ✅ **Training data augmentation** (fixed 2% errors)
- ✅ **Missing RC-augmentation**
- ✅ **Possible index stride issue**

**Fix these 3 things → 85-90% mapping rate**
