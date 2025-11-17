# GenoCache V4 - Quick Start Guide

**For:** Continuing development  
**Updated:** 2025-11-12

---

## 🚀 Project Status

- **Phase:** Week 1, Day 1 ✅ COMPLETE
- **Next:** Day 2 - Data generation
- **Goal:** 99.7% accuracy, 200-300× minimap2 speed

---

## 📁 Project Structure

```
genocache-v4/
├── training/augmentation.py  ✅ Working - NeuralAligner augmentation
├── models/encoder.py          ✅ Working - V3 multi-scale CNN
├── requirements.txt           ✅ Done
└── .venv/                     ✅ Virtual environment ready
```

---

## 💻 Quick Commands

### Activate Environment
```bash
cd /home/nebius/genocache/genocache-v4
source .venv/bin/activate
```

### Test Augmentation
```bash
python training/augmentation.py
```

### Install More Dependencies
```bash
pip install torch biopython pysam
```

---

## 📋 Day 2 Checklist

### Morning
- [ ] Download GRCh38.fa (if not exists)
- [ ] Install BioPython
- [ ] Write data generation script

### Afternoon
- [ ] Generate 130M training examples
- [ ] Validate data quality
- [ ] Save to disk (HDF5)

### Evening
- [ ] Set up train_v4.py
- [ ] Test training loop
- [ ] Start chr22 training

---

## 🎯 Key Files

| File | Purpose | Status |
|------|---------|--------|
| `training/augmentation.py` | Data augmentation | ✅ Done |
| `models/encoder.py` | Model architecture | ✅ Done |
| `training/train_v4.py` | Training loop | ⏳ Day 2 |
| `training/dataset.py` | Data loader | ⏳ Day 2 |

---

## 🔗 Important Links

- **Full Spec:** `/home/nebius/specs/2025-11-12-genocache-gpu-native-genomics-platform-full-20-week-production-plan.md`
- **V3 Lessons:** `/home/nebius/genocache/genocache-v3-backup/docs/SESSION_SUMMARY.md`
- **NeuralAligner:** `/home/nebius/genocache/genocache-v3-backup/neuraligner.md`

---

## 📊 Progress

**Week 1:** █████░░░░░░░░░░░░░░░ 20%  
**Overall:** █░░░░░░░░░░░░░░░░░░░ 1%

---

## ⚡ Key Principles

1. Build for 2027, not 2025
2. GPU-native, cloud-first
3. Validate at gates
4. Network effects are the moat

---

**Ready to continue!** 🚀
