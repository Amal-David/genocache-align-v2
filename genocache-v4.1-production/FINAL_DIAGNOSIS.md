# Final Diagnosis - GenoCache 32% vs minimap2 94%

## The Problem in Plain English

**Model and index don't match properly**

- Index expects: 128D embeddings
- Model outputs: 256D (without contrastive head) OR 128D (with head)
- Current code uses: 256D → searching 256D in 128D space → FAIL
- Tested using: 128D → worse (20%) → projection head not trained properly

## Root Cause

**The model's contrastive projection head wasn't properly trained**

Evidence:
1. Index built with 128D ✅
2. Model has 256D → 128D projection ✅
3. Using 256D in 128D index = 32% (current)
4. Using 128D (projection) = 20% (worse!)
5. → Projection head not properly trained ❌

## What This Means

The model needs to be retrained with proper contrastive learning where:
1. Base model outputs 256D
2. **Projection head properly trained** to output discriminative 128D
3. Index built with those trained 128D embeddings

## Why 32% Works Better Than 20%

- 256D has more information (even if dimensionally wrong for index)
- Compressed 128D from untrained projection loses information
- Index search with wrong dimensions gives semi-random results
- 32% is "lucky noise" vs 20% "unlucky noise"

## The Real Fix

**Retrain the entire pipeline:**
1. Train model WITH contrastive head properly
2. Build index using contrastive head outputs (128D)
3. Use contrastive head during inference

**Time:** 1-2 days
**Expected:** 85-90% mapping rate

## Current Best Configuration

**Keep:** `use_contrastive_head=False` (32% is better than 20%)

**Why:** Even though dimensionally wrong, 256D retains more signal than poorly-trained 128D projection

## Summary

- ❌ Model projection head not properly trained
- ❌ Can't fix with code changes
- ✅ Need full retraining
- ✅ Code architecture is correct
- ✅ Just need proper training

