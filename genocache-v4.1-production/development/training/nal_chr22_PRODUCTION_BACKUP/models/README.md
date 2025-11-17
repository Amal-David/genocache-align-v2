# Models Directory

The trained models are NOT included in this backup due to size (5.5 MB each).

## To Get Models

### Option 1: Train From Scratch (Recommended)
```bash
cd ../training/
python3 train_chr22_focused.py
```
This takes ~1.5 hours and creates:
- `models/chr22_nal_512bp_final.pt` (512bp seed model)
- `models/chr22_nal_256bp_final.pt` (256bp seed model, optional)

### Option 2: Copy From Original Location
```bash
cp /home/nebius/genocache/genocache-v4.1-production/development/training/nal_single_chr/models/chr22_nal_512bp_final.pt ./
```

## Model Specifications

- **Architecture:** CNN-based (4 conv layers)
- **Embedding dim:** 128
- **Seed length:** 512bp
- **Training:** 8000 batches with curriculum learning (5% → 15% errors)
- **Chromosome:** Chr22 only
- **Parameters:** ~1M

## Files Expected
- `chr22_nal_512bp_final.pt` - Main production model (REQUIRED)
- `chr22_nal_256bp_final.pt` - Alternative model (optional)
