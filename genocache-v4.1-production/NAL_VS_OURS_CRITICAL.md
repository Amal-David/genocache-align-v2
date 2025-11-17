# CRITICAL DIFFERENCE: NAL vs GenoCache

## What NeuralAligner Actually Does

### From the Paper (Line 27, 142):

**Architecture:**
```
Encoder → 128D output → Trainable Projection → InfoNCE Loss (TRAINING ONLY)
```

**Key Quote:**
> "First, for contrastive training, we remove the language modeling head and take 
> **the final encoder output as the seed representation**. To reduce the risk of 
> dimensional collapse in contrastive learning, we add **a trainable projection 
> layer that maps embeddings before they are used for the loss function**."

### What This Means:

**During TRAINING:**
- Encoder → 128D
- Projection layer → (dimension unclear, but used for loss)
- InfoNCE loss computed on projected embeddings

**During INFERENCE & INDEXING:**
- Encoder → 128D ✅
- **NO projection layer used!** ✅
- Index built with 128D encoder outputs ✅
- Search uses 128D encoder outputs ✅

**Embedding dimension: D = 128** (explicit in paper)

---

## What GenoCache Does (WRONG!)

### Current Architecture:

```
Encoder → 256D → Projection → 128D (optional)
```

**Problem:**
- Base encoder outputs **256D** (not 128D!)
- Projection head supposed to compress to 128D
- But projection not properly trained
- Index built with 128D (correct)
- Inference uses 256D (WRONG!)

---

## The ROOT CAUSE

### NAL Approach (Correct):
1. **Encoder natively outputs 128D** ✅
2. Projection only for training (dimensional collapse prevention)
3. Inference uses 128D directly from encoder
4. Index uses 128D directly from encoder
5. Everything matches!

### GenoCache Approach (Wrong):
1. **Encoder outputs 256D** ❌ (too large!)
2. Projection compresses 256D → 128D
3. Projection not properly trained
4. Index built with 128D (orphaned from training)
5. Inference uses 256D (dimension mismatch!)

---

## The REAL Problem

**We made the encoder output 256D instead of 128D!**

### Why This Happened:

**File: `genocache_core/encoder.py` line ~153**

```python
def __init__(
    self,
    emb_dim=256,  # ❌ WRONG! NAL uses 128
    seed_len=512,
    ...
):
```

**GenoCache claims:**
> "256D embeddings (vs 128D) for better specificity"

**But this breaks the design!**
- NAL: 128D native, no projection needed during inference
- GenoCache: 256D native, needs projection to 128D for index
- Projection not trained properly → 128D useless
- Forced to use 256D in 128D index → dimension mismatch!

---

## The ACTUAL Fix

### Option 1: Match NAL Architecture (BEST) ✅

**Change encoder to output 128D natively:**

```python
# File: genocache_core/encoder.py
def __init__(
    self,
    emb_dim=128,  # ← Change from 256 to 128
    seed_len=512,
    ...
):
```

**Then:**
1. Remove/ignore projection head during inference
2. Encoder directly outputs 128D
3. Matches index (128D)
4. No training issues!

### Option 2: Fix Projection Training ⚠️

**Keep 256D encoder, train projection properly:**
1. Retrain model with proper contrastive learning on 128D outputs
2. Rebuild index with trained 128D embeddings
3. Use projection head during inference

**Problem:** More complex, goes against NAL design

---

## Why NAL Works at 94%

**Simple architecture:**
- Hyena-DNA tiny: 0.5M params
- Native 128D output
- No complex projection
- Clean training: encoder → 128D → InfoNCE
- Clean inference: encoder → 128D → search
- Everything aligned!

## Why GenoCache Fails at 32%

**Over-complicated:**
- Custom encoder: 1.2M params
- 256D output (2× NAL)
- Projection layer adds complexity
- Training: encoder → 256D → projection → 128D → InfoNCE
- Inference: encoder → 256D (wrong!)
- Dimension mismatch with 128D index!

---

## The Irony

**GenoCache claims to be "Beyond NeuralAligner"**

From `encoder.py` line 3-5:
> "Improvements over NeurALigner's Hyena-tiny (0.5M params):
> - 256D embeddings (2× NeurALigner's 128D for better specificity)"

**But this "improvement" broke it!**
- Larger dimension doesn't help if mismatched
- Added complexity without benefit
- Broke the clean NAL design

---

## The Fix (Simple!)

### 1. Change Encoder Dimension (5 min)

**File:** `genocache_core/encoder.py` line ~153

```python
# BEFORE (WRONG)
def __init__(self, emb_dim=256, ...):

# AFTER (CORRECT - matches NAL)
def __init__(self, emb_dim=128, ...):
```

### 2. Retrain Model (6-8 hours)

```bash
cd development/training
python3 train_fullgenome.py --emb-dim 128 --epochs 30
```

### 3. Rebuild Index (2-3 hours)

```bash
python3 build_faiss_improved.py \
    --model models/genocache_128d.pt \
    --output indexes/genocache_128d.index
```

### 4. Update Inference (already correct!)

```python
# No change needed - encoder outputs 128D natively
embedding = self.model(x, use_contrastive_head=False)  # 128D now!
```

**Expected: 32% → 85-90% mapping rate**

---

## Summary

| Aspect | NAL (Works) | GenoCache (Broken) |
|--------|-------------|-------------------|
| Encoder output | 128D ✅ | 256D ❌ |
| Projection | Training only ✅ | Required but not trained ❌ |
| Index | 128D ✅ | 128D ✅ |
| Inference | 128D ✅ | 256D ❌ |
| Match? | YES ✅ | NO ❌ |
| Result | 94% ✅ | 32% ❌ |

**Root cause:** Encoder dimension mismatch (256D vs 128D)

**Fix:** Change emb_dim to 128, retrain, rebuild index

**Time:** 1 day

**Expected:** 85-90% mapping rate
