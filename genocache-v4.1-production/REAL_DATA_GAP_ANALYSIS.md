# Real Data Performance Gap Analysis

**Date:** 2025-11-16  
**Critical Issue:** GenoCache 21% vs minimap2 94% mapping rate on real GIAB HG002 data

---

## Test Results (100 Real GIAB HG002 Reads)

| Metric | GenoCache | minimap2 | Gap |
|--------|-----------|----------|-----|
| **Mapped** | 21% | 94% | **-73%** ❌ |
| **Unmapped** | 79% | 6% | +73% |
| **Speed** | 0.93 r/s | ~10 r/s | 11× slower |

**Critical Finding:** GenoCache fails to map 79% of real ONT reads that minimap2 successfully maps.

---

## Root Cause Analysis

### 1. Model Training Data Mismatch ⚠️ PRIMARY ISSUE

**Current State:**
- Model trained on synthetic data with simulated errors
- May not capture real ONT error patterns
- Training data distribution != real sequencing

**NeuralAligner Approach:**
- Trained with contrastive learning
- Error rates: 0.01-0.1 (uniform distribution)
- Badread simulator for realistic ONT errors
- Translation continuity enforced

**Gap:**
- Our model may lack robustness to real ONT error patterns
- Real ONT errors: homopolymer issues, systematic biases
- Synthetic training != real sequencing distribution

### 2. Index Search Parameters ⚠️ SECONDARY ISSUE

**Current:**
```python
nprobe = 64  # Number of clusters to search
top_k = 32   # Candidates per seed
```

**NeuralAligner:**
```python
nprobe = 8-32   # Fewer clusters (faster, more focused)
top_k = 32      # Same
```

**Issue:** nprobe=64 may be too high, causing:
- Retrieval of too many noisy candidates
- Slower search
- False positives diluting signal

### 3. Seed Extraction Strategy ⚠️

**Current:** Evenly spaced seeds
**NeuralAligner:** Unclear if evenly spaced or adaptive

**Potential Issue:**
- Real reads may have quality variations
- Error clustering in specific regions
- Need quality-aware seeding

### 4. Chaining Tolerance ⚠️

**Current:**
```python
colinearity_tolerance = 1000bp  # Max position error
```

**NeuralAligner:**
- Uses parameter C for positional tolerance
- Relaxed constraints for inexact anchors
- Stripe-based diagonal matching

**Issue:** May be too strict for real ONT reads with:
- Structural variations
- Large insertions/deletions
- Complex error patterns

---

## Comparison: GenoCache vs NeuralAligner

### ✅ Matches NeuralAligner

| Feature | GenoCache | NeuralAligner | Status |
|---------|-----------|---------------|--------|
| Seed count | 5 → 16 | 5 → 16 | ✅ |
| Seed length | 512bp | 256-512bp | ✅ |
| top_k | 32 | 32 | ✅ |
| Chain scoring | Count anchors | Count anchors | ✅ |
| Rescue logic | 3 conditions | Adaptive | ✅ |
| WFA alignment | Yes | Yes | ✅ |

### ⚠️ Differences from NeuralAligner

| Feature | GenoCache | NeuralAligner | Impact |
|---------|-----------|---------------|--------|
| **nprobe** | 64 | 8-32 | HIGH ⚠️ |
| **Model training** | Unknown | Contrastive + RC-aug | HIGH ⚠️ |
| **Error robustness** | Unknown | 0.01-0.1 errors | HIGH ⚠️ |
| **Translation continuity** | Unknown | Enforced | MEDIUM ⚠️ |
| **Index stride** | Unknown | ≤ seed_len/8 | MEDIUM ⚠️ |
| **Chaining algorithm** | Simple colinearity | Vectorized stripe | MEDIUM ⚠️ |

---

## Recommended Actions (Priority Order)

### 🔴 Critical (High Impact)

**1. Reduce nprobe to 16-32** (Quick Fix - 5 min)
```python
# In adaptive_seeding.py
if hasattr(self.index, 'nprobe'):
    self.index.nprobe = 16  # Start conservative (NAL uses 8-32)
```
**Expected:** 20-30% mapping rate improvement

**2. Relax Chaining Tolerance** (Quick Fix - 10 min)
```python
# Increase from 1000bp to 3000-5000bp
colinearity_tolerance: int = 3000  # More lenient for real reads
```
**Expected:** 10-20% mapping rate improvement

**3. Lower Minimum Chain Score** (Quick Fix - 5 min)
```python
# In rescue logic
if chains[0].score < self.min_seeds / 3:  # Was /2
    need_rescue = True
```
**Expected:** 5-10% improvement

### 🟡 Important (Medium Impact)

**4. Verify Model Training** (Investigation - 30 min)
- Check how model was trained
- Verify error rates used
- Compare with NAL training protocol
**Expected:** Understanding of root cause

**5. Analyze Failed Reads** (Investigation - 1 hour)
- Why did 79% fail to map?
- Embedding quality issues?
- Search returning no candidates?
- Chain formation failing?
**Expected:** Targeted fixes

### 🟢 Optional (Long-term)

**6. Retrain Model** (If needed - 1-2 days)
- Use NAL training protocol
- Contrastive learning with 0.01-0.1 errors
- RC-augmentation
- Translation continuity

**7. Implement NAL Chaining** (Enhancement - 3-4 hours)
- Vectorized stripe-based chaining
- Better handling of inexact anchors
- Equation 2 from paper

---

## Quick Validation Plan (2 hours)

### Phase 1: Parameter Tuning (30 min)

1. **Test nprobe values:**
   - Try: 16, 24, 32
   - Run on same 100 GIAB reads
   - Measure mapping rate

2. **Test colinearity_tolerance:**
   - Try: 2000, 3000, 5000
   - Run on same reads
   - Measure mapping rate

3. **Test rescue threshold:**
   - Try: min_seeds/3, min_seeds/4
   - Run on same reads
   - Measure mapping rate

### Phase 2: Analyze Failures (1 hour)

1. **Instrument code to track:**
   - How many seeds get candidates?
   - How many chains form?
   - Where do most reads fail?

2. **Check embedding quality:**
   - Are embeddings discriminative?
   - Do similar sequences cluster?
   - Distance distribution analysis

### Phase 3: Compare Best Config (30 min)

- Run best configuration on 100 reads
- Compare with minimap2
- Document improvement

**Target:** 50-70% mapping rate (from 21%)

---

## Hypothesis Testing

### Hypothesis 1: nprobe Too High
**Test:** Run with nprobe=16 vs 64  
**Expected:** Higher mapping rate, better precision  
**Time:** 10 minutes

### Hypothesis 2: Chaining Too Strict
**Test:** Increase tolerance from 1000 to 5000bp  
**Expected:** More chains form, higher mapping rate  
**Time:** 10 minutes

### Hypothesis 3: Model Not Robust
**Test:** Check embedding distances for failed reads  
**Expected:** Large distances, poor discrimination  
**Time:** 30 minutes  
**Solution:** Model retraining required

---

## Expected Outcomes

### Conservative (Parameter tuning only)
- **Mapping rate:** 21% → 40-50%
- **Gap to minimap2:** 94% - 50% = 44%
- **Status:** Better but not production-ready

### Optimistic (Parameters + minor code fixes)
- **Mapping rate:** 21% → 60-70%
- **Gap to minimap2:** 94% - 70% = 24%
- **Status:** Acceptable for research use

### Ideal (Model retraining required)
- **Mapping rate:** 21% → 85-90%
- **Gap to minimap2:** 94% - 90% = 4%
- **Status:** Production-ready

---

## Decision Matrix

| Action | Time | Expected Gain | Risk |
|--------|------|---------------|------|
| **Reduce nprobe to 16** | 5 min | +20-30% | LOW |
| **Increase tolerance** | 10 min | +10-20% | LOW |
| **Lower rescue threshold** | 5 min | +5-10% | LOW |
| **Analyze failures** | 1 hour | Understanding | NONE |
| **Code fixes** | 2-3 hours | +10-20% | MEDIUM |
| **Model retraining** | 1-2 days | +30-60% | HIGH |

**Recommended Path:**
1. Quick parameter tuning (20 min) → 40-50% mapping rate
2. Analyze failures (1 hour) → Understand root cause
3. Decide: Accept 50% OR invest in retraining for 85%+

---

## Summary

**Current:** 21% mapping rate on real data ❌  
**Issue:** Likely model training + parameter mismatch  
**Quick fixes:** Can reach 40-50% in 20 minutes  
**Proper fix:** Model retraining → 85%+ (1-2 days)  

**Recommendation:** Try quick fixes first, then decide on retraining investment.
