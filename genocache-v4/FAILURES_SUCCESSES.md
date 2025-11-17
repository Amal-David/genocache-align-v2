# GenoCache V4 - Failures & Successes Log

**Purpose:** Track what works and what doesn't for rapid iteration

---

## ✅ Successes

### Day 1 (2025-11-12)
1. **Augmentation pipeline** - Implemented cleanly, all tests pass
2. **V3 architecture reuse** - Perfect fit, no changes needed
3. **Project structure** - Clean organization

---

## ❌ Failures & Solutions

### Day 1 (2025-11-12)
1. **Missing numpy** 
   - Error: ModuleNotFoundError
   - Solution: Created venv, installed dependencies ✅

2. **Genome loading failed (0 chromosomes)** 
   - Error: ValueError: empty range in randrange(0, -511)
   - Root cause: Filter expected `chr1` format but GRCh38 uses `NC_000001.11`
   - Solution: Updated filter to accept both `NC_` and `chr` prefixes ✅
   - Result: Successfully loaded 25 chromosomes, 3.08 billion bp ✅

---

3. **Variable sequence lengths in DataLoader**
   - Error: RuntimeError: Trying to resize storage that is not resizable
   - Root cause: Augmentation created sequences of different lengths
   - Solution: Force all sequences to exactly 512bp (pad or trim) ✅
   - Result: Training started successfully ✅

---

## ⚠️ Current Blockers

**None** - Training running successfully!

---

**Last Updated:** 2025-11-12
