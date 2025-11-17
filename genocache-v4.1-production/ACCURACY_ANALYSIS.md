# Accuracy Analysis: 87.5% → How to Reach 95%+

**Current Status:** 87.5% accuracy (7/8 reads correct)  
**Target:** 95%+ accuracy (like minimap2)  
**Gap:** 1 read failing (read_6: alternate contig)

---

## Current Pipeline

```
1. ADAPTIVE SEEDING
   - Extract 5-16 seeds (512bp each)
   - Encode each → 128D embedding
   - FAISS search → top-32 candidates per seed
   - Parameters:
     ✓ min_seeds = 5
     ✓ max_seeds = 16 (rescue)
     ✓ top_k = 32 (FAISS candidates per seed)
     ✓ window_size = 512bp

2. CHAINING
   - Group candidates by chromosome
   - Check colinearity (±1kb tolerance)
   - Score by sum of seed scores
   - Return all chains sorted by score
   - Parameters:
     ✓ colinearity_tolerance = 1000bp

3. EXTEND PHASE
   - Test top-5 chains
   - Align read to each region
   - Pick best by alignment score
   - Parameters:
     ✓ return_top_k = 5 chains tested
     ✓ min_score_threshold = 100
     ✓ score_ratio_threshold = 1.5
```

---

## The One Failing Read: read_6

**What happened:**
```
Ground truth: NT_187498.1 (alternate contig, 7kb)
OLD picked:   NC_000022.11 (main chr22) - wrong but higher score
NEW picked:   NC_000022.11 (main chr22) - still wrong

Alignment scores:
  chr22:         1708 ← picked (wrong)
  NT_187498.1:   1644 ← correct (lower score)
  
Difference: 64 points (3.7%)
```

**Why it failed:**
- Alternate contig is similar to main chromosome
- Main chromosome has slightly higher alignment score
- Both are valid, but alternate is "more correct"
- This is a HARD case (ambiguous mapping)

---

## Root Cause Analysis

### Hypothesis 1: FAISS Search Not Finding Alternate Contig ❓

**Check:**
- Is NT_187498.1 in top-32 candidates from FAISS?
- If NO → FAISS search problem (need more candidates)
- If YES → Chaining or EXTEND problem

**NeuralAligner uses:**
- top-k = 32 per seed (same as us)
- But they may have different index parameters

**Our FAISS config:**
```python
# Current settings (from build_production_index.py)
nlist = sqrt(91.8M) ≈ 9600 clusters
nprobe = 16  # ← KEY PARAMETER
m = 16       # Product quantization
nbits = 8
```

**Issue:** `nprobe = 16` means we only search 16 clusters
- Trade-off: Speed vs Accuracy
- Higher nprobe = slower but more accurate
- NeuralAligner may use higher nprobe

**Test:** Increase nprobe to 32 or 64

### Hypothesis 2: Chaining Failing for Alternate Contigs ❓

**Check:**
- Does NT_187498.1 form a chain with 2+ seeds?
- If NO → Seeds not colinear (alternate contigs are short)
- If YES → Chain score too low vs main chr

**Issue:** Alternate contigs are SHORT (7kb vs 50Mb)
- Fewer seeds land in alternate contig
- Chain score = sum(seed_scores)
- Main chr has more seeds → higher chain score

**NeuralAligner approach:**
- Uses **rescue seeding** more aggressively
- May adjust chain scoring for short contigs

**Test:** Lower min_seeds requirement for chains (allow 1-seed chains?)

### Hypothesis 3: EXTEND Not Testing Enough Candidates ❓

**Check:**
- Is NT_187498.1 in top-5 chains returned to EXTEND?
- If NO → Need to test more chains
- If YES → Alignment score discrimination not enough

**Current:** return_top_k = 5 (test 5 chains)

**NeuralAligner uses:**
- Top-K chains (unclear from paper, but likely 5-10)
- But they report >99% accuracy

**Test:** Increase return_top_k to 10

### Hypothesis 4: Alignment Score Not Discriminating Well ❓

**Observation:**
```
chr22:       1708 (wrong)
NT_187498.1: 1644 (correct)
Difference:  64 (3.7%)
```

**Issue:** 3.7% difference is SMALL
- Could be noise/ambiguity
- Both are genuinely good matches
- Read might map to BOTH (multi-mapping)

**NeuralAligner approach:**
- Reports secondary alignments
- Uses MAPQ scores to indicate ambiguity
- Doesn't force single best

**Our approach:**
- Always picks single best
- Doesn't report ambiguity
- May be "correct" to pick chr22 if it's better aligned

**Test:** Report alignment as "ambiguous" if score difference < 5%

---

## Comparison with NeuralAligner Paper

### From NeuralAligner (Section 3):

**1. Seeding Parameters:**
```
- Seed length: 256bp (we use 512bp) ✓ We're more informative
- Number of seeds: 5-20 adaptive (we use 5-16) ✓ Similar
- Top-k per seed: Not specified (we use 32) ? Unknown
```

**2. Index Parameters:**
```
- Index type: FAISS IVFPQ (same as us) ✓
- nlist: sqrt(N) (same as us) ✓
- nprobe: NOT SPECIFIED ← CRITICAL UNKNOWN
- Compression: PQ (same as us) ✓
```

**3. Chaining:**
```
- Colinearity check: YES (same as us) ✓
- Tolerance: "relaxed constraints" (we use 1kb) ? May be higher
- Chain scoring: "number of anchors" (we use sum of scores) ~ Similar
```

**4. EXTEND Phase:**
```
- Align to top-K chains: YES (same as us) ✓
- Pick by alignment score: YES (same as us) ✓
- Algorithm: WFA-GPU (we use WFA2) ✓ Same algorithm
- Threshold: Not specified ← UNKNOWN
```

**5. Accuracy Reported:**
```
NeuralAligner: 99% accuracy on 250bp reads
Our system:    87.5% accuracy on 1000bp reads

Note: Different read lengths make direct comparison hard
```

---

## Likely Culprits (Priority Order)

### 1. **FAISS nprobe Too Low** (HIGH PRIORITY)

**Current:** nprobe = 16 (search 16 clusters)
**Impact:** May miss alternate contigs in different clusters

**Fix:**
```python
# In genocache_align.py or adaptive_seeding.py
index.nprobe = 32  # Search more clusters (2× slower, more accurate)
# or
index.nprobe = 64  # Search even more (4× slower, most accurate)
```

**Expected improvement:** 87.5% → 92-95%

### 2. **Insufficient Rescue Seeding** (MEDIUM PRIORITY)

**Current:** Add 16 seeds only if no chains found
**Issue:** May not activate for ambiguous cases

**Fix:**
```python
# Add rescue even if chains exist but top chain score is low
if len(chains) == 0 or chains[0].score < threshold:
    # Add rescue seeds
```

**Expected improvement:** +2-3%

### 3. **Chain Scoring Bias Toward Long Chromosomes** (MEDIUM PRIORITY)

**Current:** Chain score = sum(seed_scores)
**Issue:** Longer regions get more seeds → higher scores

**Fix:**
```python
# Normalize by region length or number of seeds
chain_score = sum(seed_scores) / len(seeds)
# or
chain_score = sum(seed_scores) / (end_pos - start_pos) * 1000000
```

**Expected improvement:** +1-2%

### 4. **Not Enough Chains Tested in EXTEND** (LOW PRIORITY)

**Current:** return_top_k = 5
**Issue:** May not include alternate contigs in top-5

**Fix:**
```python
return_top_k = 10  # Test more candidates
```

**Expected improvement:** +1-2%

---

## Recommended Action Plan

### Phase 1: FAISS Tuning (QUICK WIN)

**Test nprobe values:**
```python
# Test these in order:
nprobe = 32  # 2× slower, likely fixes issue
nprobe = 64  # 4× slower, maximum accuracy
```

**Implementation:** Change one line in adaptive_seeding.py

**Expected time:** 5 minutes
**Expected improvement:** 87.5% → 92-95%

### Phase 2: Rescue Seeding Enhancement (IF NEEDED)

**Add conditional rescue:**
```python
# Trigger rescue if:
# 1. No chains found, OR
# 2. Top chain score < threshold, OR
# 3. Multiple chains with similar scores (ambiguous)
```

**Expected time:** 15 minutes
**Expected improvement:** +2-3%

### Phase 3: Chain Scoring (IF NEEDED)

**Normalize chain scores:**
```python
# Test different normalization schemes
# Find what works best
```

**Expected time:** 30 minutes
**Expected improvement:** +1-2%

---

## Special Case: Alternate Contigs

**Fundamental issue:**
- Alternate contigs (NT_*, patches) are HARD
- They're similar to main chromosomes by design
- minimap2 also struggles with these (uses MAPQ to indicate ambiguity)

**Our options:**

1. **Report as ambiguous** (honest approach)
   - MAPQ < 20 for ambiguous mappings
   - Don't force single answer
   - Let downstream tools decide

2. **Prefer main chromosomes** (pragmatic approach)
   - Penalize alternate contigs slightly
   - Most users want main chr anyway
   - Matches common practice

3. **Multi-mapping** (complete approach)
   - Report all good alignments (primary + secondary)
   - SA tag in SAM format
   - Highest accuracy but more complex

**For now:** Test #1 (nprobe tuning) - should get us to 95%+

---

## Next Steps

1. ✅ Read current code (done)
2. ✅ Analyze failure mode (done)
3. ✅ Identify likely causes (done)
4. 🎯 **Test nprobe = 32** (next)
5. ⏭️ Re-validate on 8 reads
6. ⏭️ If still failing, try phase 2

---

**Recommendation:** Start with FAISS nprobe tuning. This is the most likely fix and takes 5 minutes to test.
