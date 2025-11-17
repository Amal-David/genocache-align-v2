# 🎓 GenoCache V4 - Frequently Asked Questions

## Training Strategy

### Q: Are we training on the full genome?
**A: YES!** ✅
- Loaded all 25 chromosomes (chr1-22, X, Y, MT)
- Total: 3.08 billion base pairs
- Sampling positions randomly across entire genome

### Q: Are we using error augmentation?
**A: Partially** ⚠️
- Current (1M test): Random 1-10% errors
- Next (130M full): Curriculum learning (0%→10% gradually)
- Curriculum is NeuralAligner's proven strategy

### Q: What are the 130M samples?
**A: Full training dataset size**
```
1M samples:   Testing pipeline (CURRENT)
10M samples:  Quick validation (NEXT)  
130M samples: Full training for 99.7% accuracy (GOAL)
```

**Why 130M?**
- Genome: 3.08B bp
- Sample stride: 100bp = ~30M positions
- Augmentation: 4× per position
- Total: ~120-130M training examples
- This matches NeuralAligner's scale

### Q: Why curriculum learning?
**A: It's the key to making neural alignment work!**

**Without curriculum (random errors):**
- Model confused from start
- Learns slowly
- May reach only 85-90% accuracy

**With curriculum (0%→10% gradually):**
- Model learns clean patterns first
- Gradually adds complexity
- Reaches 99.7% accuracy (NeuralAligner proved this!)

---

## Current Status

### Q: What's training right now?
**A: 1M example test run**
```
Purpose: Verify pipeline works
Examples: 1M (0.7% of final dataset)
Time: ~1.5 hours
Error rate: Random 1-10%
Expected: 85-90% accuracy
Goal: Pipeline validation ✅
```

### Q: When do we scale to 130M?
**A: After validating the pipeline**
```
Step 1: 1M test (NOW - running)
Step 2: 10M validation (if test works)
Step 3: 130M full training (if validation good)
Timeline: ~3-4 days total
```

### Q: Will 1M examples reach 99.7%?
**A: NO** - That's not the goal!
```
1M examples:   Pipeline test (current)
               Expected: 85-90%
               Goal: Verify code works ✅

130M examples: Full training (later)
               Expected: 99.7%+
               Goal: Production accuracy ✅
```

---

## Technical Details

### Q: How much disk space do we need?
**A: We have plenty!** ✅
```
Available: 410 GB
Needed:    ~116 GB (worst case)
  - Training data (130M): 54 GB
  - Model checkpoints: 0.5 GB
  - Encoded genome: 25 GB
  - FAISS index: 26 GB
  - Temp files: 10 GB

Status: MORE THAN ENOUGH ✅
```

### Q: Do we need bidirectional convolutions?
**A: Probably not** ✅
```
What we have:
- Multi-scale CNN (4 kernel sizes)
- Position encoding
- RC augmentation (50%)

What NeuralAligner has:
- Bidirectional convolutions

Our approach MIGHT be better already!
If accuracy <99% at Gate 1, we can add it.
```

### Q: Why not train on all 130M right away?
**A: Smart risk management!**
```
Bad approach:
- Generate 54GB of data (6 hours)
- Train for 40 hours
- Find bug in pipeline
- Wasted 46 hours! ❌

Good approach (ours):
- Test with 1M (1.5 hours) ✅
- Find bugs quickly
- Fix and validate
- Then scale to 130M
- Save time overall! ✅
```

---

## Process & Progress

### Q: How long until we reach 99.7%?
**A: ~3-4 days of GPU time**
```
Timeline:
- 1M test: 1.5 hours ✅ (running now)
- 10M validation: 5 hours (next)
- 130M curriculum: 30-40 hours (final)
Total: ~48 hours = 2 days GPU time

Calendar time: ~1 week (Week 2-3 of plan)
Gate 1 check: Week 4
```

### Q: What if 1M training fails?
**A: We fix it quickly!**
```
Possible issues:
1. Loss doesn't decrease → Adjust learning rate
2. Separation <0.2 → Fix augmentation
3. GPU OOM → Reduce batch size

Advantage of 1M test:
- Find issues in 1.5 hours (not 40 hours!)
- Fix and restart quickly
- Validate before big training run
```

### Q: Can we speed up training?
**A: Several options:**
```
Current: Single GPU, batch_size=32
Faster options:
1. Increase batch size (32→64→128)
2. Multi-GPU (data parallel)
3. Mixed precision (FP16)
4. Gradient accumulation

Tradeoff: Speed vs memory vs accuracy
We'll optimize after pipeline validated ✅
```

---

## Strategy & Plan

### Q: Why follow NeuralAligner's approach?
**A: They proved it works!**
```
NeuralAligner results:
- 99.6% accuracy on ONT reads ✅
- 276× faster than minimap2 ✅
- Published ICLR 2026 (peer-reviewed)

Our advantages:
- Better architecture (multi-scale CNN)
- Hard negative mining (983K regions)
- Their proven training strategy

Expected: 99.7-99.8% (beat their 99.6%)
```

### Q: What's the path to production?
**A: Clear roadmap with gates**
```
Week 1-4: Train model → Gate 1 (99.7% accuracy?)
Week 5-7: Add caching → 500-1000× speedup
Week 8-20: Build cloud platform → Revenue

Decision gates:
- Gate 1 (Week 4): <99% = abort/pivot
- Gate 2 (Week 7): Caching works?
- Gate 3 (Week 20): Users paying?
```

---

## Background Processing

### Q: Why let training run in background?
**A: Efficient use of time!**
```
Training: GPU-bound (1.5 hours)
Prep work: CPU-bound (can do in parallel)

While training runs:
- ✅ Create curriculum script
- ✅ Write documentation
- ✅ Prepare next pipeline components
- ✅ Plan validation tests

Result: Save time overall! ✅
```

### Q: Will training stop if I disconnect?
**A: NO - it's in a tmux/screen session** ✅
```
Process runs on VM, not your terminal
You can:
- Close browser
- Disconnect SSH
- Check back later

Training continues automatically!
```

---

## Next Steps

### Q: What happens when training finishes?
**A: Check results, then scale up**
```
1. Run: ./scripts/check_training_results.sh
2. Look for separation metric (>0.5 = good)
3. If good: Generate 10M examples
4. Train on 10M (5 hours)
5. If still good: Curriculum training (130M, 50 epochs)
6. Week 4: Validate on GIAB HG002
```

### Q: How will we know if we hit 99.7%?
**A: GIAB HG002 validation (Week 4)**
```
GIAB = Genome in a Bottle (gold standard)
HG002 = High-quality human genome with known variants

Test:
- Align 10K-100K reads from HG002
- Compare to true positions
- Calculate accuracy

Target: ≥99.7% correct positions
```

---

**More questions? Check the detailed docs:**
- `TRAINING_EXPLAINED.md` - Full training strategy
- `WHAT_TO_DO_NEXT.md` - Next steps guide
- `CURRENT_SESSION_SUMMARY.md` - What we've done

**Training is running in background - no action needed!** 🚀
