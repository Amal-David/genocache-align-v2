# FINAL VALIDATION - Complete Instructions

**Purpose:** Validate EXTEND phase fixes 37% → 95%+ accuracy  
**Data:** All test files located and ready  
**Status:** ⏳ Awaiting PyTorch environment

---

## COMPLETE TEST PLAN

### Test 1: Reproduce Bug (Already Done)
**Result:** 37.5% chromosome accuracy confirmed ✅

```bash
# Already have these files:
test_complete_10reads.sam      # OLD method (37% accuracy)
minimap2_same_10reads.sam      # Baseline (100% accuracy)

# Comparison already shows bug:
python3 compare_chromosome_accuracy.py
# Output: 37.5% chromosome accuracy (3/8 correct)
```

### Test 2: Validate Mock Fix (Already Done)
**Result:** 0% → 100% improvement ✅

```bash
python3 test_extend_mock.py
# Output: OLD 0%, NEW 100%, +100% improvement
```

### Test 3: Real Data with EXTEND Fix (TO DO)
**Goal:** Prove 37% → 95%+ on actual reads

```bash
cd /home/nebius/genocache/genocache-v4

# Step 1: Run NEW pipeline with EXTEND phase
python3 run_with_extend.py \
    --input test_10_reads_exact.fa \
    --output test_extend_10reads.sam \
    --model models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt \
    --index indexes/genocache_v4_production.index \
    --metadata indexes/genocache_v4_production.metadata.pkl

# Step 2: Run minimap2 for comparison (already have this)
# minimap2_same_10reads.sam already exists

# Step 3: Compare ALL three
python3 compare_three_way.py \
    --old test_complete_10reads.sam \
    --new test_extend_10reads.sam \
    --minimap minimap2_same_10reads.sam \
    --output FINAL_COMPARISON_REPORT.md

# Expected output:
# OLD (seed count):  37.5% accuracy
# NEW (EXTEND):      95%+ accuracy  ← PROOF OF FIX!
# minimap2:          100% accuracy
```

### Test 4: Large Scale Validation
**Goal:** Validate on 500 reads

```bash
# Run on larger dataset
python3 run_with_extend.py \
    --input /home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa \
    --output test_extend_500reads.sam \
    --model models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt

# Compare with minimap2
minimap2 -ax map-ont \
    /home/nebius/genocache/GRCh38.fa \
    /home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa \
    > minimap2_500reads.sam

# Analyze results
python3 analyze_accuracy.py \
    --genocache test_extend_500reads.sam \
    --minimap minimap2_500reads.sam \
    --output LARGE_SCALE_RESULTS.md
```

---

## SCRIPTS TO CREATE

### Script 1: run_with_extend.py

```python
#!/usr/bin/env python3
"""
Run complete pipeline with EXTEND phase
"""

import torch
import faiss
import pickle
from Bio import SeqIO
from adaptive_seeding import AdaptiveSeeder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner
from models.encoder import GenoCacheEncoder

def main():
    # Load model
    model = GenoCacheEncoder(emb_dim=128, seed_len=512)
    checkpoint = torch.load('models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt')
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    # Load index
    index = faiss.read_index('indexes/genocache_v4_production.index')
    with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
        metadata = pickle.load(f)
    
    # Initialize components
    seeder = AdaptiveSeeder(model, index, metadata)
    aligner = FastAligner()
    extender = ExtendPhase(aligner)
    
    # Process reads
    with open('test_extend_10reads.sam', 'w') as out:
        for record in SeqIO.parse('test_10_reads_exact.fa', 'fasta'):
            # Step 1: Get top-k candidates
            candidates = seeder.align_read(str(record.seq), record.id, return_top_k=5)
            
            if candidates:
                # Step 2: EXTEND - align to each, pick best by score
                result = extender.extend_and_score(str(record.seq), candidates)
                
                if result:
                    # Write SAM entry
                    write_sam_entry(out, record.id, result)

if __name__ == '__main__':
    main()
```

### Script 2: compare_three_way.py

```python
#!/usr/bin/env python3
"""
Three-way comparison: OLD vs NEW vs minimap2
"""

def compare_three_way(old_sam, new_sam, minimap_sam):
    # Parse all three SAM files
    old_aligns = parse_sam(old_sam)
    new_aligns = parse_sam(new_sam)
    mm2_aligns = parse_sam(minimap_sam)
    
    # Ground truth from minimap2
    common_reads = set(old_aligns.keys()) & set(new_aligns.keys()) & set(mm2_aligns.keys())
    
    old_correct = 0
    new_correct = 0
    
    print("=" * 80)
    print("THREE-WAY COMPARISON")
    print("=" * 80)
    
    for read_id in sorted(common_reads):
        old_chr = old_aligns[read_id]['chr']
        new_chr = new_aligns[read_id]['chr']
        truth_chr = mm2_aligns[read_id]['chr']
        
        old_match = (old_chr == truth_chr)
        new_match = (new_chr == truth_chr)
        
        if old_match:
            old_correct += 1
        if new_match:
            new_correct += 1
        
        status = ""
        if not old_match and new_match:
            status = "✅ FIXED!"
        elif old_match and new_match:
            status = "✅ Both correct"
        elif not old_match and not new_match:
            status = "❌ Both wrong"
        elif old_match and not new_match:
            status = "⚠️ Regression"
        
        print(f"{read_id}:")
        print(f"  Ground truth: {truth_chr}")
        print(f"  OLD (seed ct): {old_chr} {'✓' if old_match else '✗'}")
        print(f"  NEW (EXTEND):  {new_chr} {'✓' if new_match else '✗'}")
        print(f"  {status}")
        print()
    
    total = len(common_reads)
    
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"OLD method (seed count): {old_correct}/{total} ({old_correct/total*100:.1f}%)")
    print(f"NEW method (EXTEND):     {new_correct}/{total} ({new_correct/total*100:.1f}%)")
    print(f"Improvement: {new_correct - old_correct} reads (+{(new_correct-old_correct)/total*100:.1f}%)")
    
    if new_correct > old_correct:
        print()
        print("✅ SUCCESS: EXTEND phase improves accuracy!")
    else:
        print()
        print("❌ FAILURE: No improvement")

if __name__ == '__main__':
    compare_three_way(
        'test_complete_10reads.sam',
        'test_extend_10reads.sam',
        'minimap2_same_10reads.sam'
    )
```

---

## EXPECTED RESULTS

### Before EXTEND Fix (Current):
```
Test: 10 reads from chr22
Method: Seed count discrimination
Results:
  - Chromosome accuracy: 37.5% (3/8 correct)
  - Position accuracy: ~980bp (when chr correct)
  - Speed: 6.7 reads/sec
Status: ❌ BUG - Too inaccurate
```

### After EXTEND Fix (Expected):
```
Test: Same 10 reads
Method: Alignment score discrimination (EXTEND phase)
Results:
  - Chromosome accuracy: 95%+ (8/8+)  ← FIXED!
  - Position accuracy: ~980bp (unchanged)
  - Speed: 1-3 reads/sec (slower but accurate)
Status: ✅ WORKING - Production ready
```

### Large Scale (500 reads):
```
Test: 500 synthetic chr22 reads
Method: EXTEND phase
Results:
  - Chromosome accuracy: 90-95%
  - Position accuracy: <1kb median
  - Speed: ~500 sec total
  - Errors: Mostly repeat regions
Status: ✅ Validated on large dataset
```

---

## SUCCESS CRITERIA

### Minimum (Proof):
- [ ] NEW method > 70% accuracy (vs OLD 37%)
- [ ] Clear improvement demonstrated
- [ ] No regressions on correct reads

### Target (Production):
- [ ] NEW method > 90% accuracy
- [ ] Comparable to minimap2 (within 10%)
- [ ] Position accuracy maintained

### Ideal:
- [ ] NEW method > 95% accuracy
- [ ] Validated on 500+ reads
- [ ] Error analysis complete

---

## TIMELINE

### Quick Test (1 hour):
1. Fix PyTorch environment (30 min)
2. Run on 10 reads (15 min)
3. Compare results (15 min)

### Full Validation (4 hours):
1. Quick test (1 hour)
2. Run on 500 reads (2 hours)
3. Generate comprehensive report (1 hour)

### With Report (5 hours):
1. Full validation (4 hours)
2. Write up results (1 hour)
3. Create presentation (optional)

---

## FILES NEEDED

### Existing:
- ✅ `test_10_reads_exact.fa` - 10 test reads
- ✅ `test_complete_10reads.sam` - OLD results (37%)
- ✅ `minimap2_same_10reads.sam` - Baseline
- ✅ `reads_chr22_synth_1kb_500.fa` - 500 reads
- ✅ `extend_phase.py` - EXTEND implementation
- ✅ `adaptive_seeding.py` - Modified for top-k
- ✅ `fast_alignment.py` - Aligner

### To Create:
- ⏳ `run_with_extend.py` - NEW pipeline script
- ⏳ `compare_three_way.py` - Comparison script
- ⏳ `analyze_accuracy.py` - Statistics script

---

## BLOCKERS

1. **PyTorch Environment**
   - Status: Not available in system Python
   - Impact: Can't load model
   - Solution: Install or use pre-extracted embeddings

2. **Model Loading**
   - File exists: `models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt`
   - Size: ~500MB
   - Needs: PyTorch to load

3. **Index Access**
   - File exists: `indexes/genocache_v4_production.index`
   - Size: ~7GB
   - Needs: faiss-cpu package (likely installed)

---

## WORKAROUND (No PyTorch)

Can validate EXTEND logic without full pipeline:

```python
# Use pre-computed candidates from OLD run
# Just test EXTEND phase scoring logic
# Compare OLD picking (seed count) vs NEW (align score)
# This proves the concept without needing the model
```

---

## NEXT SESSION COMMANDS

```bash
# 1. Check environment
cd /home/nebius/genocache/genocache-v4
python3 -c "import torch; print('PyTorch OK')" || echo "Need PyTorch"

# 2. If PyTorch available, run validation
python3 run_with_extend.py

# 3. Compare results
python3 compare_three_way.py

# 4. Check improvement
# Expected: 37% → 95%+ accuracy

# 5. Generate report
python3 generate_final_report.py
```

---

**Status:** Complete plan ready, awaiting execution  
**Confidence:** High (mock test already proved concept)  
**Expected Duration:** 1-5 hours depending on scope  
**Expected Result:** 37% → 95%+ accuracy improvement ✅
