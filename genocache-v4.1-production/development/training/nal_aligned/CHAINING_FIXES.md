# Chaining Algorithm Fixes - Following NAL Algorithm 1

**Date:** 2025-11-16  
**Issue:** Chaining implementation deviated from NAL paper Algorithm 1

---

## What Was Wrong

### Old Implementation (Dict Clustering)

```python
# OLD: Brute-force dict clustering
candidates = {}
for anchor in anchors:
    expected_y = y_ij - x_i  # Correct formula
    # Scan all dict keys to find match
    for y in list(candidates.keys()):
        if abs(y - expected_y) <= tolerance:
            candidates[y].append(anchor)
            break
```

**Problems:**
1. ❌ No sorting by Y
2. ❌ No similarity D usage
3. ❌ Dict key scanning (inefficient)
4. ❌ Wrong complexity O(#anchors × #clusters)

---

## What Was Fixed

### New Implementation (NAL Algorithm 1)

```python
# NEW: Follow Algorithm 1 exactly
# Step 1: Adjust Y for strand
for anchor in anchors:
    if strand == '+':
        adjusted_y = y_ij - x_i
    else:
        adjusted_y = y_ij + x_i + seed_len - read_len

# Step 2: Sort by adjusted Y
adjusted_anchors.sort(key=lambda a: a['adjusted_y'])

# Step 3: Merge neighbors within C_b
for i in range(len(adjusted_anchors) - 1):
    j = i + 1
    if adjusted_anchors[j]['adjusted_y'] - adjusted_anchors[i]['adjusted_y'] <= C_b:
        # Use max similarity as stripe axis
        if adjusted_anchors[j]['similarity'] > adjusted_anchors[i]['similarity']:
            adjusted_anchors[i]['adjusted_y'] = adjusted_anchors[j]['adjusted_y']
            adjusted_anchors[i]['similarity'] = adjusted_anchors[j]['similarity']
        
        # Accumulate count
        adjusted_anchors[i]['score'] += adjusted_anchors[j]['score']
        
        # Track all merged anchors
        adjusted_anchors[i]['merged_anchors'].extend(adjusted_anchors[j]['merged_anchors'])
        
        # Mark j as inactive
        active[j] = False

# Step 4: Collect active chains
# Step 5: Sort by score descending
```

**Fixes:**
1. ✅ Sort anchors by adjusted Y
2. ✅ Merge neighbors in genomic order
3. ✅ Use similarity D to pick stripe axis
4. ✅ Track all merged anchors
5. ✅ Correct complexity O(NMK log(MK))

---

## Key Changes

### 1. Sorting (Algorithm 1, line 5)

**Before:** No sorting, dict insertion order  
**After:** `adjusted_anchors.sort(key=lambda a: a['adjusted_y'])`

### 2. Similarity Usage (Algorithm 1, line 8)

**Before:** Ignored similarity scores  
**After:** 
```python
if adjusted_anchors[j]['similarity'] > adjusted_anchors[i]['similarity']:
    # Use j's position and similarity as stripe axis
    adjusted_anchors[i]['adjusted_y'] = y_j
    adjusted_anchors[i]['similarity'] = adjusted_anchors[j]['similarity']
```

### 3. Sequential Merging (Algorithm 1, lines 7-12)

**Before:** Scan all dict keys (O(n²))  
**After:** Merge sorted neighbors (O(n))

### 4. Anchor Tracking

**Before:** Lost merged anchors  
**After:** Track all anchors in each stripe
```python
adjusted_anchors[i]['merged_anchors'].extend(adjusted_anchors[j]['merged_anchors'])
```

---

## What's Still Correct

✅ **Expected Y formulas** (Equation 2):
- Positive: `y = y_ij - x_i`
- Negative: `y = y_ij + x_i + seed_len - read_len`

✅ **Score by anchor count** (not span)

✅ **Separate strands** (+/-)

✅ **Top-K selection**

---

## Parameter Documentation

### Tolerance (C and C_b)

**Paper distinction:**
- C: Colinearity tolerance in Eq. 2
- C_b: Bias tolerance for merging in Alg. 1

**Our implementation:**
- Use same value for both (reasonable simplification)
- Default: 1000bp (appropriate for long reads)

### Rescue Thresholds

**Paper:** Qualitative description ("low score", "ambiguous")

**Our design:**
- Low score: `best_score < num_seeds/2`
- Ambiguous: `second_score >= best_score × 0.8`
- Perfect: `best_score == num_seeds`

These are reasonable heuristics not explicitly specified in paper.

---

## Expected Impact

### Before Fix (71% mapping)

Issues:
- Inefficient dict clustering
- No similarity weighting
- Order-dependent behavior

### After Fix (Expected: 75-80%)

Improvements:
- Proper sorted merging
- Best-similarity stripe selection
- Deterministic O(NMK log(MK))
- More accurate chain identification

**Combined with WFA:** Expected 85-90% (matches paper)

---

## Testing

Run validation to verify improvement:

```bash
cd /home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned

source /home/nebius/genocache/.venv/bin/activate

python3 align_nal.py \
    --model models/genocache_nal.pt \
    --index indexes/genocache_nal_stride32.index \
    --positions indexes/genocache_nal_stride32.positions.npz \
    --reference /home/nebius/genocache/GRCh38.fa \
    --reads /home/nebius/genocache/genocache-v4.1-production/validation/data/giab_hg002_100reads.fastq \
    --output results/giab_nal_fixed.sam \
    --device cuda
```

Expected: 75-80% mapping (vs 71% before fix)

---

## Summary

**Fixed:** Chaining now follows NAL Algorithm 1 exactly
- ✅ Sort by adjusted Y
- ✅ Use similarity for stripe axis
- ✅ Sequential neighbor merging
- ✅ Proper O(NMK log(MK)) complexity

**Documented:** Design choices for thresholds
- Tolerance: C = C_b (simplification)
- Rescue: Explicit thresholds (reasonable)

**Ready:** For validation testing with improved accuracy
