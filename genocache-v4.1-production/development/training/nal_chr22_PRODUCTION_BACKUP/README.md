# Chr22 NAL Training - Production Backup

**Date:** 2025-11-16
**Status:** WORKING - 76.5% accuracy (beats minimap2 at 75.5%)

## What's In This Folder

This is a CLEAN backup of the working Chr22 NAL training pipeline that achieved
76.5% accuracy on synthetic reads with ground truth.

### Directory Structure

```
nal_chr22_PRODUCTION_BACKUP/
├── training/           # Training scripts
│   ├── train_chr22_focused.py       # Main training (8000 batches, curriculum)
│   └── encoder_nal.py               # Model architecture
├── indexing/           # Index building
│   └── build_chr22_index.py         # Build FAISS index for Chr22
├── alignment/          # Alignment pipeline
│   ├── seeding_chr22.py             # Seed extraction & FAISS search
│   ├── chaining_nal.py              # Anchor chaining
│   └── test_chr22_nal.py            # Complete test pipeline
├── models/             # Trained models (not included - too large)
│   └── README.md                    # Instructions to download/train
├── results/            # Test results
│   └── chr22_test.log               # Final test results (76.5%)
└── documentation/      # All docs
    ├── TRAINING_DIFF.md             # OLD vs NEW changes
    ├── RESULTS.md                   # Final results summary
    └── S3_BACKUP.md                 # S3 backup instructions

## Key Changes From Previous Version

1. **8000 batches** (vs ~100 before) - 80x MORE training
2. **Curriculum learning** - Progressive 5% → 15% error rate
3. **Chr22 only** - Focused training on single chromosome
4. **Result:** 43% → 76.5% accuracy improvement!

## How to Use

### Training
```bash
cd training/
python3 train_chr22_focused.py
# Takes ~1.5 hours, creates models/chr22_nal_512bp_final.pt
```

### Indexing
```bash
cd indexing/
python3 build_chr22_index.py
# Takes ~10 minutes, creates indexes/chr22_nal_512bp_stride32.{index,positions.npz}
```

### Testing
```bash
cd alignment/
python3 test_chr22_nal.py
# Tests on 200 synthetic Chr22 reads
```

## Results

- **Accuracy:** 76.5% (vs minimap2's 75.5%)
- **Improvement:** +33.5% vs previous full-genome model (43%)
- **Speed:** ~35ms/read (2-tier adaptive strategy)

## Requirements

- Python 3.12+
- PyTorch
- FAISS-GPU
- NumPy

