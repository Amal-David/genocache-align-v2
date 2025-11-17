# CRITICAL: NAL Accuracy Issues Detected

## 🚨 **URGENT: False Positive Rate 88.89%**

Date: 2025-11-16  
Discovery: Comparison with minimap2 on 100 GIAB HG002 reads  
Status: ❌ **PRODUCTION BLOCKED** - Major accuracy issues detected

---

## Summary

When comparing NAL alignments with minimap2 (gold standard) on the same 100 reads:

| Metric | Value | Assessment |
|--------|-------|------------|
| **Precision** | **11.11%** | ❌ **CRITICAL** |
| **False Positive Rate** | **88.89%** | ❌ **CRITICAL** |
| Recall | 91.67% | ✅ Good |
| Mapping Rate (NAL) | 99% | ✅ High |
| Mapping Rate (minimap2) | 94.1% | Reference |

**Finding:** NAL maps 99% of reads, but only 11% of those mappings are correct!

---

## Details

### Confusion Matrix

```
                    minimap2 Mapped    minimap2 Unmapped
NAL Mapped                 11                88  (FALSE POSITIVES)
NAL Unmapped                1                 0
```

### Types of Errors

1. **Wrong Chromosome:** 77 reads (87.5% of false positives)
   - NAL maps to completely different chromosome than minimap2
   - Example: NAL→Chr2, minimap2→Chr18

2. **Wrong Position:** 82 reads (93.2% of false positives)
   - Position differences >1000bp (some >100Mbp!)
   - Example: Off by 166 million base pairs

3. **Wrong Strand:** Many reads
   - NAL and minimap2 disagree on strand orientation

---

## Example Case Study

**Read:** `b095fff0-6e9d-4e30-92e8-d10ef3eac27e`

| Source | Chromosome | Position | MAPQ | Strand |
|--------|------------|----------|------|--------|
| NAL | NC_000002.12 (Chr2) | 239,776,577 | 30 | Forward |
| minimap2 | NC_000018.10 (Chr18) | 73,838,260 | 60 | Reverse |

**Discrepancy:** Completely different chromosomes, positions 166M bp apart!

---

## Root Cause Analysis

### Potential Issues

1. **Index Building Problem** ⚠️ HIGH LIKELIHOOD
   - Seeds may be indexed to wrong genomic positions
   - Position calculation error during index creation
   - Stride or offset bug in build_index_nal.py

2. **Reference Mismatch** ⚠️ MEDIUM LIKELIHOOD
   - Training used different reference than alignment
   - Reference version inconsistency (GRCh38 vs hg38)

3. **Chaining Algorithm Bug** ⚠️ MEDIUM LIKELIHOOD
   - Picking incorrect chains from anchor set
   - Strand orientation error
   - Position calculation error

4. **Seed Position Calculation** ⚠️ HIGH LIKELIHOOD
   - Error in converting FAISS index positions back to genomic coordinates
   - Off-by-one errors or stride miscalculation
   - Chromosome boundary handling

5. **WFA Alignment Artifacts** ⚠️ LOW LIKELIHOOD
   - WFA is just final alignment, doesn't affect initial mapping
   - But could affect reported position

---

## Evidence

### What Works
- ✅ Training converged (loss 0.19)
- ✅ Index built successfully (81.6M seeds)
- ✅ Seeds retrieved from FAISS (99% efficiency)
- ✅ Chains formed successfully
- ✅ WFA generates proper CIGAR strings
- ✅ High mapping rate (99%)

### What's Broken
- ❌ Mapped positions are wrong (88.89% incorrect)
- ❌ Wrong chromosomes selected
- ❌ Position errors >100Mbp in some cases
- ❌ Strand orientation errors

---

## Why High Mapping Rate Misled Us

**The Trap:**
- NAL reports **99-100% mapping rate**
- All internal metrics looked good (chain scores, anchor density)
- We assumed high mapping = correct mapping
- **But:** High mapping rate just means "found a chain", not "correct location"

**The Reality:**
- NAL is confidently mapping reads to wrong locations
- High MAPQ scores (15-60) suggest confidence
- But minimap2 disagrees 89% of the time
- **This is worse than low mapping - it's wrong answers!**

---

## Debugging Steps Required

### Immediate (Critical)

1. **Verify Index Building**
   ```bash
   # Check if positions.npz has correct genomic coordinates
   python3 -c "
   import numpy as np
   pos_data = np.load('indexes/genocache_nal_stride32.positions.npz', allow_pickle=True)
   print('Chr names:', pos_data['chr_names'][:10])
   print('Positions shape:', pos_data['positions'].shape)
   print('Sample positions:', pos_data['positions'][:10])
   "
   ```

2. **Trace Single Read Through Pipeline**
   - Extract one failing read
   - Print seeds extracted
   - Print FAISS results (indices + positions)
   - Print converted genomic positions
   - Compare with minimap2 result
   - Find where it goes wrong

3. **Check Reference Consistency**
   ```bash
   # Verify training and alignment use same reference
   md5sum GRCh38.fa  # Used for alignment
   # Compare with reference used during training
   ```

4. **Validate Position Calculation**
   ```python
   # In seeding_nal.py, verify:
   anchor['ref_chr'] = self.chr_names[chr_idx]
   anchor['ref_pos'] = positions[i]  # Is this correct?
   ```

### Secondary (Important)

5. **Test on Known-Good Read**
   - Align read that minimap2 maps to Chr1:1000000
   - Check where NAL maps it
   - If also wrong, confirms systematic issue

6. **Check for Off-By-One Errors**
   - Stride calculation
   - Position indexing (0-based vs 1-based)
   - Chromosome boundaries

7. **Validate Against Paper**
   - Re-check all position calculations vs NAL paper
   - Verify our interpretation matches theirs

---

## Impact Assessment

### Current State
❌ **NAL cannot be used for production**  
❌ **All previous results are unreliable**  
❌ **99-100% mapping claims are misleading**

### What This Means

1. **Previous Test Results Invalid**
   - 99-100% mapping was measuring "found any chain"
   - Not "found correct location"
   - Parameter sweeps still valid (relative comparisons)
   - But absolute accuracy claims are wrong

2. **Need to Fix Before Proceeding**
   - Cannot do 1000-read test until fixed
   - Cannot deploy adaptive tiers until fixed
   - Cannot claim superior accuracy until fixed

3. **Comparison with minimap2**
   - NAL: 99% mapped, 11% correct = **11% effective accuracy**
   - minimap2: 94% mapped, ~90%+ correct = **~85%+ effective accuracy**
   - **minimap2 is 7-8x more accurate!**

---

## Recovery Plan

### Phase 1: Root Cause (1-2 days)
1. Trace failing read through pipeline
2. Identify exact point of failure
3. Fix position calculation bug
4. Verify fix on 10 reads

### Phase 2: Validation (1 day)
1. Re-run 100-read test
2. Compare with minimap2
3. Target: <5% false positive rate
4. Target: >90% precision

### Phase 3: Full Testing (2-3 days)
1. Test on 1000 reads
2. Compare with minimap2
3. Validate adaptive tiers
4. Document fixes

---

## Lessons Learned

### What Went Wrong

1. **Trusted Internal Metrics Too Much**
   - Chain scores looked good
   - Anchor density looked good
   - But didn't compare with ground truth early enough

2. **Didn't Validate Early**
   - Should have compared with minimap2 at 71% mapping stage
   - Could have caught this much earlier
   - Wasted time optimizing wrong system

3. **Misunderstood Success Metrics**
   - High mapping rate ≠ correct mapping
   - Need precision/recall, not just recall
   - Should always compare with gold standard

### What To Do Differently

1. **Always Compare with Gold Standard**
   - Compare with minimap2 from day 1
   - Don't trust internal metrics alone
   - Validate on every major change

2. **Start Simple, Validate Early**
   - Test one read first
   - Verify it's correct before scaling
   - Don't optimize until baseline works

3. **Check Assumptions**
   - Verify position calculations match paper
   - Test edge cases (chromosome boundaries, etc.)
   - Don't assume "it worked" means "it's correct"

---

## Next Steps

1. **STOP all optimization work** ❌
2. **Debug position calculation** ⚠️
3. **Fix root cause** 🔧
4. **Re-validate everything** ✅
5. **Only then proceed with tiers/optimization**

---

## Status

**BLOCKED:** Cannot proceed until accuracy issues resolved.

**Action Required:**
- Debug position calculation in index building or anchor retrieval
- Verify all coordinate transformations
- Compare single read trace with minimap2
- Fix and re-test

**Timeline:** Estimated 3-5 days to fix and validate.

---

## Files for Investigation

Priority order:
1. `build_index_nal.py` - Index building (position storage)
2. `seeding_nal.py` - Position retrieval from index
3. `chaining_nal.py` - Position reporting in chains
4. `align_nal.py` - Final position in SAM output

Check each for:
- Position calculation errors
- Stride/offset bugs
- Coordinate system mismatches
- Off-by-one errors

---

**Report Date:** 2025-11-16  
**Severity:** CRITICAL  
**Priority:** P0 - Block all other work  
**Owner:** Requires immediate attention
