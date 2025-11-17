# Autonomous Validation Pipeline - Implementation Plan

**Status:** Running autonomously  
**Goal:** Fix 37% chromosome accuracy bug and validate on real data  
**Time:** ~6-8 hours autonomous execution

---

## Current Status

✅ **EXTEND phase implemented** - Modified adaptive_seeding.py to return top-k candidates  
✅ **Score-based discrimination implemented** - extend_phase.py created  
✅ **WFA-GPU built** - libwfagpu.so compiled successfully  
✅ **Bug confirmed** - 37.5% chromosome accuracy measured  
✅ **Root cause identified** - Picking by seed count instead of alignment score

⚠️ **Blockers:**
- PyTorch not available in system Python (need to work around)
- Need to complete WFA-GPU Python bindings
- Need to download GIAB HG002 data

---

## Autonomous Execution Plan

### Phase 1: Complete WFA-GPU Integration (1-2 hours)

**Step 1.1: Create complete Python wrapper**
- Use ctypes to wrap WFA-GPU C API
- Implement batch alignment interface
- Test on sample sequences

**Step 1.2: Integrate with FastAligner**
- Modify fast_alignment.py to use WFA-GPU
- Fallback to parasail if WFA-GPU unavailable
- Benchmark speed improvement

### Phase 2: Test Fixed Pipeline (30 min)

**Step 2.1: Create minimal test (no PyTorch needed)**
- Use pre-computed candidates (mock seeding)
- Test EXTEND phase with WFA-GPU
- Verify chromosome accuracy improvement

**Step 2.2: Compare with minimap2**
- Run on same 10 reads
- Measure chromosome accuracy (expect 37% → 95%+)
- Measure position accuracy
- Measure speed

### Phase 3: Download GIAB HG002 Data (30 min)

**Step 3.1: Download ONT reads**
```bash
# GIAB HG002 ONT data (chr22 subset)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_ucsc_ONT_chr22.fastq.gz
```

**Step 3.2: Extract subset**
- Take first 100 reads
- Length distribution: 1-100kb
- Representative of real ONT data

### Phase 4: Full Validation (2-3 hours)

**Step 4.1: Run GenoCache on GIAB data**
- Full pipeline with EXTEND phase
- Generate SAM output
- Measure speed (reads/sec)

**Step 4.2: Run minimap2 on same reads**
- Standard minimap2 -ax map-ont
- Generate SAM output for comparison
- Measure speed

**Step 4.3: Compare results**
- Chromosome-level accuracy
- Position-level accuracy (median error)
- CIGAR agreement
- Speed comparison

### Phase 5: Generate Report (30 min)

**Step 5.1: Create comprehensive comparison**
- Accuracy metrics (before/after EXTEND)
- Speed metrics (GenoCache vs minimap2)
- Error analysis (where do we still fail?)
- Production readiness assessment

**Step 5.2: Document results**
- Update HOW_NEURALIGNER_SOLVED_IT.md
- Create VALIDATION_RESULTS.md
- Create summary for hackathon

---

## Alternative Approaches (If Blockers Occur)

### If PyTorch unavailable:
1. **Option A:** Use pre-extracted embeddings from index
   - Load vectors directly from FAISS
   - Skip model loading
   - Still test EXTEND phase

2. **Option B:** Mock seeding results
   - Create synthetic candidates
   - Test EXTEND phase logic only
   - Validate score-based discrimination

3. **Option C:** Install PyTorch in user space
   ```bash
   python3 -m pip install --user torch --index-url https://download.pytorch.org/whl/cu118
   ```

### If WFA-GPU bindings difficult:
1. **Use parasail for now**
   - 250× slower but works
   - Proves EXTEND phase concept
   - Migrate to WFA-GPU later

2. **Direct subprocess call**
   - Use WFA-GPU command-line tool
   - Parse output
   - Less elegant but functional

### If GIAB download fails:
1. **Use synthetic reads**
   - Generate from chr22 reference
   - Known ground truth
   - Still validates pipeline

2. **Use existing test data**
   - Expand current 10-read set to 100
   - Less realistic but proves concept

---

## Success Criteria

### Minimum (Proof of Concept):
- ✅ EXTEND phase implemented
- ✅ Chromosome accuracy improves (37% → 70%+)
- ✅ Tested on 10+ reads
- ✅ WFA-GPU or parasail working

### Target (Production Ready):
- ✅ Chromosome accuracy 90%+ on GIAB data
- ✅ Position accuracy <5kb median error
- ✅ Speed: 1-5 reads/sec with EXTEND
- ✅ WFA-GPU integrated for production

### Stretch (Publication Quality):
- ✅ Chromosome accuracy 95%+ on GIAB
- ✅ Position accuracy <1kb median
- ✅ Comprehensive comparison with minimap2
- ✅ Error analysis and failure modes documented

---

## Timeline

**Hour 0-2:** WFA-GPU bindings + integration  
**Hour 2-3:** Test fixed pipeline, compare with minimap2  
**Hour 3-4:** Download GIAB data  
**Hour 4-7:** Full validation on GIAB  
**Hour 7-8:** Generate reports and documentation

**Total:** 8 hours autonomous execution

---

## Current Task: WFA-GPU Python Bindings

**Status:** In progress  
**Next:** Complete ctypes wrapper, test on sample sequences

**Files being created:**
- `wfa_gpu_wrapper.py` - Python bindings
- `test_wfagpu_alignment.py` - Test script
- `autonomous_validation.py` - Full pipeline

---

**Last Updated:** 2025-11-14 07:00 UTC  
**Mode:** AUTONOMOUS - Continuing without intervention
