# 📋 What to Do When Training Completes

**Training ETA:** ~1.5 hours from start (around 23:45 UTC)

---

## ✅ Step 1: Check Training Results

```bash
cd /home/nebius/genocache/genocache-v4

# View final results
tail -100 logs/training_initial.log

# Look for these lines:
# - "Val positive similarity: 0.XXX"
# - "Val negative similarity: 0.XXX"
# - "Separation: 0.XXX"
```

### ✅ Success Criteria:
- Positive similarity: >0.8
- Negative similarity: <0.3
- Separation: >0.5

---

## 🎯 Step 2A: If Training Succeeded (likely)

### Scale up training data:
```bash
# Generate 10M examples (~4.2GB)
.venv/bin/python training/generate_training_data.py \
  --num-examples 10000000 \
  --batch-size 100000 \
  --output data/training_10M.h5

# Train on 10M
.venv/bin/python training/train_simple.py \
  --data data/training_10M.h5 \
  --epochs 10 \
  --checkpoint models/checkpoints/model_10M.pt
```

### Then implement curriculum learning:
- Error rate progression: 0% → 1% → 3% → 5% → 8% → 10%
- 50 epochs total
- This is the NeuralAligner strategy

---

## 🔧 Step 2B: If Training Failed (unlikely)

### Debug checklist:
1. **Check loss didn't diverge** (should decrease, not increase)
2. **Check separation** (if <0.2, model didn't learn)
3. **Possible fixes:**
   - Lower learning rate (1e-4 → 1e-5)
   - Increase batch size (32 → 64)
   - Add more regularization (dropout)
   - Check data quality (some examples may be corrupted)

---

## 📈 Step 3: Continue to Week 4 (Gate 1)

### Timeline:
- **Day 3:** Train on 10M examples
- **Day 4-5:** Implement curriculum learning (50 epochs)
- **Day 6-7:** Build alignment pipeline + test on chr22
- **Week 4:** Full GIAB validation → Check if 99.7% accuracy

### Gate 1 Decision:
- ✅ 99.7% accuracy → Proceed to Phase 2 (caching)
- ⚠️ 95-99% → Add improvements (curriculum, bidirectional)
- ❌ <95% → Major debugging needed

---

## 💾 Storage Status

**Current usage:** 5.6 GB  
**Available:** 410 GB  
**Enough for:** Full 130M training dataset ✅

---

## 🚀 Current Status

```
Week 1: ████████░░░░░░░░░░ 40% (AHEAD OF SCHEDULE)

Day 1: ████████████████████ 100% ✅ Complete
Day 2: ████████████░░░░░░░░  60% 🔥 Training...
Day 3: Pending
Week 4: Gate 1 validation
```

---

**Don't forget to check `logs/training_initial.log` when training completes!**
