# Validation Plan - Real Dataset Testing

**Date:** 2025-11-15  
**Purpose:** Validate improvements on real sequencing data (not synthetic)

---

## Dataset Selection

Based on NeuralAligner paper and available public data:

### Option 1: GIAB HG002 ONT Reads (RECOMMENDED)

**Source:** Genome in a Bottle Consortium  
**Location:** `s3://ont-open-data/giab_2025.01/`  
**Latest Release:** 2025.01  

**Specifications:**
- Sample: HG002 (PGP Ashkenazi Son)
- Technology: Oxford Nanopore (PromethION)
- Basecalling: Latest model (Dorado)
- Quality: High-quality, production-grade data

**What we'll use:**
- Chromosome 22 reads only (manageable size)
- ~100-1000 reads for validation
- Mix of read lengths (typical ONT distribution)

**Download command:**
```bash
# Using AWS CLI (if available)
aws s3 cp s3://ont-open-data/giab_2025.01/HG002/ . --recursive --no-sign-request

# Or using wget/curl with public URLs
```

### Option 2: Badread-Simulated (Backup)

If real data download fails:
```bash
# Generate using Badread (NeuralAligner's method)
badread simulate \
    --reference GRCh38_chr22.fa \
    --quantity 100x \
    --length 5000,15000 \
    --identity 90,99,100 \
    --error_model nanopore2023 \
    > simulated_chr22_reads.fastq
```

---

## Validation Protocol

### Phase 1: Baseline (Original Code)

1. Backup current modified code
2. Restore original `adaptive_seeding.py` (before improvements)
3. Run on real dataset (100 reads)
4. Record accuracy

### Phase 2: With Improvements

1. Restore modified code (chain scoring + rescue logic)
2. Run on same dataset
3. Record accuracy
4. Compare with baseline

### Phase 3: Analysis

1. Identify which reads improved
2. Categorize errors (main chr vs alternates, etc.)
3. Statistical significance test
4. Document findings

---

## Directory Structure

```
genocache-v4.1-production/
├── genocache_core/
│   ├── adaptive_seeding_BASELINE.py     # Original (before changes)
│   ├── adaptive_seeding_WITH_IMPROVEMENTS.py  # Modified version
│   └── adaptive_seeding.py              # Current (symlink or copy)
│
├── validation/
│   ├── data/
│   │   ├── giab_hg002_chr22_100reads.fastq
│   │   └── ground_truth_chr22.txt
│   ├── results/
│   │   ├── baseline_results.txt
│   │   └── improved_results.txt
│   └── scripts/
│       ├── download_giab_data.sh
│       ├── run_baseline.py
│       └── run_improved.py
│
└── docs/
    ├── CHANGES_LOG.md                   # All code changes
    ├── VALIDATION_PLAN.md               # This file
    └── VALIDATION_RESULTS.md            # Results summary
```

---

## Metrics to Track

### Primary Metrics

1. **Accuracy:** Percentage of correctly mapped reads
2. **Precision:** True positives / (True positives + False positives)
3. **Recall:** True positives / (True positives + False negatives)

### Secondary Metrics

4. **Speed:** Reads per second
5. **Alternate contig handling:** Specifically track performance on alternates
6. **Ambiguous cases:** Count of low-confidence mappings

### Error Analysis

7. **Error types:**
   - Wrong chromosome (main)
   - Wrong chromosome (alternate)
   - Unmapped (should map)
   - Mapped (should be unmapped)

---

## Expected Outcomes

### Baseline (Original Code)

- Accuracy: 85-90% (estimated)
- Issues: Biased toward main chromosomes

### With Improvements

- Accuracy: 88-93% (estimated +3-5%)
- Better: Alternate contig handling
- Better: Ambiguous case resolution

### Statistical Test

- Use paired t-test or McNemar's test
- Significance level: p < 0.05
- Report confidence intervals

---

## Timeline

1. **Download data:** 30 min - 1 hour
2. **Prepare baseline:** 15 min
3. **Run baseline validation:** 1-2 hours
4. **Run improved validation:** 1-2 hours
5. **Analysis & report:** 1 hour

**Total:** 4-6 hours

---

## Success Criteria

✅ **Success if:**
- Improvement > 2% accuracy (statistically significant)
- No regression on main chromosomes
- Better alternate contig handling

⚠️ **Investigate if:**
- No improvement (<1%)
- Regression on any metric
- High variance in results

❌ **Rollback if:**
- Accuracy decreases
- Major bugs discovered
- Performance significantly worse

---

## Documentation Requirements

After validation, document:

1. **Exact dataset used**
   - Version, source, size
   - Number of reads
   - Read length distribution

2. **Results**
   - Baseline vs improved (tables)
   - Statistical significance
   - Error analysis

3. **Conclusions**
   - What improved?
   - What didn't?
   - Next steps

---

## Next Steps

1. ⏭️ Create validation directory structure
2. ⏭️ Download GIAB HG002 chr22 data
3. ⏭️ Create baseline copy of original code
4. ⏭️ Run validation (baseline + improved)
5. ⏭️ Analyze and document results

---

**Status:** READY TO START  
**Priority:** HIGH  
**Estimated Time:** 4-6 hours
