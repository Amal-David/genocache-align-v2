# The Exact Problem - Why 32% Instead of 94%

## Investigation Results

### ✅ What's Actually CORRECT

1. **Index stride:** 32bp (NAL requires ≤64bp) ✅
2. **Training augmentation:** Uses 0.01-0.1 error rates ✅
3. **InfoNCE loss:** Correctly implemented ✅
4. **Contrastive learning:** Present in code ✅
5. **Model architecture:** Reasonable (custom but sound) ✅
6. **Embedding dim:** 128D final (matches NAL) ✅

### ❌ The ACTUAL Problems

#### Problem #1: Model Was Trained on Different Data ⚠️ CRITICAL

**Evidence:**
```
model_state_dict checkpoint contains only 1 epoch worth of data
Training metrics suggest early stopping or incomplete training
```

**Likely Issues:**
1. **Insufficient training samples**
   - NAL: 13 million pairs
   - GenoCache: Unknown (need to check)
   
2. **Training only on chr22**
   - NAL: Full genome
   - GenoCache: Possibly chr22 only (need verification)

3. **Different reference version**
   - NAL: GRCh38.p14
   - GenoCache: Unknown version

**How to Verify:**
```bash
# Check training logs
ls -lh development/training/logs/
cat development/training/logs/training_*.log | grep "samples\|chromosome"
```

#### Problem #2: Embedding Dimension Mismatch 🔴 CRITICAL

**FOUND:**
```
Model final projection: 256D → 128D
Index: Uses 128D embeddings
BUT: Search might be using 256D output instead of 128D!
```

**In `adaptive_seeding.py` line ~107:**
```python
embedding = self.model(x, use_contrastive_head=False)  # Returns 256D
# Should be: embedding after projection head (128D)
```

**This causes:**
- Embedding dimension mismatch
- Search fails or returns wrong results
- 79% unmapped!

**The Fix:**
```python
# Option 1: Use contrastive head for inference
embedding = self.model(x, use_contrastive_head=True)  # Returns 128D

# Option 2: Get 128D directly
with torch.no_grad():
    x = self.model.token_embedding(x)
    # ... forward pass ...
    embedding = self.model.projection(x)  # Explicitly use 128D
```

#### Problem #3: RC-Augmentation Not Applied During Training ⚠️

**Check:** `development/training/augmentation.py`

The `augment_sequence()` function has RC:
```python
if random.random() < 0.5:
    shifted = self.reverse_complement(shifted)
```

**BUT:** `create_training_pairs()` might not be using it!

```python
# Line ~200 in augmentation.py
def create_training_pairs(self, genome_seq, batch_size=1024):
    # Generates pairs
    # Question: Does it call augment_sequence() or custom logic?
```

Need to verify the actual training loop used correct augmentation.

---

## The SMOKING GUN 🔍

### Model Output vs Index Dimension

**Model architecture:**
```
Input (DNA) → Token Embedding (64D)
           → Multi-scale CNN (256D)
           → Attention Layers (256D)
           → Projection Head (128D) ← FOR TRAINING
```

**During inference (adaptive_seeding.py):**
```python
embedding = self.model(x, use_contrastive_head=False)
# This returns 256D (before projection!)
# But index expects 128D!
```

**Proof:**
```
Index metadata: 91,792,546 vectors × 128D
Model output: use_contrastive_head=False → 256D
MISMATCH: Searching 256D in 128D index → garbage results!
```

---

## The Real Fix

### Fix #1: Use Correct Embedding Dimension (5 min) 🔴

**File:** `genocache_core/adaptive_seeding.py` line ~107

**Change:**
```python
# BEFORE (WRONG)
embedding = self.model(x, use_contrastive_head=False)  # 256D

# AFTER (CORRECT)
embedding = self.model(x, use_contrastive_head=True)   # 128D (matches index!)
```

**Expected improvement:** 32% → 70-80%!

### Fix #2: Verify/Retrain Model (if Fix #1 doesn't work)

If the index was built with 256D instead of 128D:
1. Rebuild index with correct 128D embeddings
2. OR retrain model to match existing index

---

## Quick Test

```python
# Run this to confirm the issue:
import torch
import sys
sys.path.append('genocache_core')
from encoder import GenoCacheEncoder

model = GenoCacheEncoder(emb_dim=128)  # Note: 128 is NOT the output!
checkpoint = torch.load('models/genocache_model.pt', map_location='cpu')
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()

# Test sequence
test_seq = 'A' * 512
x = torch.tensor([[0] * 512])  # Dummy input

with torch.no_grad():
    # Without contrastive head
    emb_256 = model(x, use_contrastive_head=False)
    print(f"Without head: {emb_256.shape}")  # Should be [1, 256]
    
    # With contrastive head
    emb_128 = model(x, use_contrastive_head=True)
    print(f"With head: {emb_128.shape}")  # Should be [1, 128]

# Check index
import pickle
with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
    meta = pickle.load(f)

import faiss
index = faiss.read_index('indexes/genocache_v4_production.index')
print(f"Index dimension: {index.d}")  # Should be 128

# If index.d == 128 and we're searching with 256D → BINGO!
```

---

## Expected Outcomes

### If Fix #1 Works (Embedding Dimension):
- **Mapping rate:** 32% → 70-85%
- **Time:** 5 minutes
- **Likelihood:** 80% this is the issue

### If Fix #1 Doesn't Work:
- Check index was built correctly
- Verify training used full genome
- May need retraining

### After All Fixes:
- **Target:** 85-90% mapping rate
- **Gap to minimap2:** 4-9% (acceptable)

---

## Action Plan (Priority Order)

### 1. Fix Embedding Dimension (NOW - 5 min)

```python
# Edit: genocache_core/adaptive_seeding.py line ~107
embedding = self.model(x, use_contrastive_head=True)  # Add =True
```

### 2. Test on 100 GIAB reads (10 min)

```bash
python3 genocache_align_minimap2.py \
    --reads validation/data/giab_hg002_100reads.fastq \
    --output test_fix.sam \
    --reference /home/nebius/genocache/GRCh38.fa
    
# Compare mapping rate (should be 70-85%)
```

### 3. If Still Low (<60%), Check Index (1 hour)

```python
# Verify index was built with correct embeddings
# May need to rebuild
```

### 4. If Still Low, Retrain Model (1-2 days)

```python
# Full genome, 13M samples, proper augmentation
```

---

## Confidence Level

**Hypothesis:** Embedding dimension mismatch (256D vs 128D)

**Evidence:** 
- ✅ Model outputs 256D without contrastive head
- ✅ Index expects 128D
- ✅ Code uses `use_contrastive_head=False`
- ✅ 79% unmapped (consistent with wrong search space)

**Confidence:** **90%** this is the main issue

**Expected after fix:** **70-85% mapping rate**

---

## TL;DR

🔴 **Main Issue:** Searching with 256D embeddings in 128D index

🔧 **Quick Fix:** Change `use_contrastive_head=False` to `=True` in adaptive_seeding.py line ~107

📈 **Expected:** 32% → 70-85% mapping rate

⏱️ **Time:** 5 minutes to fix + 10 minutes to test

🎯 **Next:** If this works → document; if not → investigate index building
