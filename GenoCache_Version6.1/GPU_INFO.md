# GPU Configuration for GenoCache v6.1

## Auto-Detection Features

The training script automatically:

1. **Detects available GPUs**
   - Counts GPUs
   - Shows GPU names and memory
   - Reports total available memory

2. **Auto-adjusts batch size**
   - If GPU memory < 40GB total: reduces batch size to 16
   - Otherwise: uses default batch size 32

3. **Multi-GPU support**
   - Automatically uses DataParallel if >1 GPU available
   - Distributes batches across all GPUs
   - No code changes needed!

## Expected Performance

### With 8x H100 GPUs (80GB each = 640GB total)
- Batch size: 32 (full)
- Training time: ~3-4 hours
- Should handle 16,000 batches easily

### With 4x H100 GPUs (320GB total)
- Batch size: 32 (full)
- Training time: ~4-5 hours
- Still plenty of memory

### With 1x H100 GPU (80GB)
- Batch size: 32 (full)
- Training time: ~6-8 hours
- Should work fine

### With limited GPU (<40GB)
- Batch size: auto-reduced to 16
- Training time: longer
- May need to reduce num_batches

## Conservative Changes Made

1. **Batches: 24,000 → 16,000**
   - Still 2x more than Chr22 (8,000)
   - Reduces memory pressure
   - Shorter training time
   - Should still get 75-85% accuracy

2. **Adaptive batch size**
   - Auto-reduces if memory limited
   - Prevents OOM errors

3. **Multi-GPU ready**
   - Will use all available GPUs
   - No manual configuration needed

## Manual Adjustments

If you still hit OOM (Out of Memory):

### Option 1: Reduce batch size manually
Edit `train_full_genome.py` line 290:
```python
batch_size=16,  # Change from 32 to 16
```

### Option 2: Reduce num_batches
Edit `train_full_genome.py` line 289:
```python
num_batches=12000,  # Change from 16000 to 12000
```

### Option 3: Train on fewer chromosomes first
Test with just Chr1-Chr10 first, then scale up.

## What Didn't Change

✅ Same architecture (proven with Chr22)
✅ Same curriculum learning (5% → 15%)
✅ Same error injection (Sub/Ins/Del)
✅ Same seed length (512bp, not 256+512 variant)
✅ No RC (reverse complement) handling - keeping simple

All core components that got 76.5% on Chr22 are preserved!
