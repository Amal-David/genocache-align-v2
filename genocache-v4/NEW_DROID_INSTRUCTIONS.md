# Instructions for New Droid Session - 8-GPU Training

**Context:** Training GenoCache V4 with optimal batch size (8192) using 8 GPUs  
**Goal:** Achieve 99%+ accuracy (vs current 96.6%)  
**Time:** ~5-8 hours for full training

---

## Quick Start for New Droid

Hi! You're taking over 8-GPU training for GenoCache V4. Here's what you need to know:

### Background
- **Current status:** Model trained with batch=1024 (single GPU, memory limited)
- **Current accuracy:** 96.6% @ ±1kb tolerance
- **Goal:** Train with batch=8192 (8 GPUs) to reach 99%+ accuracy
- **Why:** InfoNCE loss works best with large batches (more negative samples)

### What to Do

**Step 1: Download everything from S3**
```bash
# User will provide S3 bucket name
export S3_BUCKET="user-provided-bucket-name"

# Download the download script
aws s3 cp s3://${S3_BUCKET}/genocache-v4/scripts/download_from_s3.sh .
chmod +x download_from_s3.sh

# Edit to set bucket name
sed -i "s/YOUR_BUCKET_NAME/${S3_BUCKET}/g" download_from_s3.sh

# Run it
./download_from_s3.sh
```

**Step 2: Setup environment**
```bash
# Create conda environment
conda create -n genocache python=3.11 -y
conda activate genocache

# Install PyTorch with CUDA
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install dependencies
pip install biopython==1.81 numpy==1.24.3 tqdm==4.65.0 \
    faiss-gpu==1.7.4 pandas==2.0.3 matplotlib==3.7.2 seaborn==0.12.2

# Verify GPU access
python3 -c "import torch; print(f'GPUs: {torch.cuda.device_count()}')"
# Should print: GPUs: 8
```

**Step 3: Read the setup guide**
```bash
cat SETUP_GUIDE.md
# This has complete instructions, but here's the TL;DR
```

**Step 4: Launch training**
```bash
# Test on single GPU first (5 minutes)
python3 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 1 \
  --output-dir ./test

# If that works, launch full 8-GPU training
nohup torchrun --nproc_per_node=8 scripts/train_8gpu.py \
  --genome GRCh38.fa \
  --batch-size 1024 \
  --epochs 50 \
  --lr 1e-4 \
  --output-dir ./models/8gpu_batch8192 \
  > training.log 2>&1 &

# Get PID
echo $! > training.pid
```

**Step 5: Monitor training**
```bash
# Watch logs
tail -f training.log

# Check GPU usage
watch -n 5 nvidia-smi

# Check if still running
ps -p $(cat training.pid)
```

**Step 6: Upload results back to S3**
```bash
# During training (sync every hour)
watch -n 3600 'aws s3 sync ./models/ s3://${S3_BUCKET}/genocache-v4-8gpu/models/'

# After training completes
aws s3 sync ./models/ s3://${S3_BUCKET}/genocache-v4-8gpu/models/
aws s3 cp training.log s3://${S3_BUCKET}/genocache-v4-8gpu/logs/
```

---

## What to Expect

**Training Progress:**
- **Per epoch:** ~6-10 minutes (50 epochs total)
- **Total time:** ~5-8 hours
- **Key metric:** Separation (pos_sim - neg_sim)
  - Current best: 12.40
  - Target: 13-15+

**Success Criteria:**
- Separation > 13.0 (better than current)
- Positive similarity > 0.95
- Negative similarity < 0.00
- Loss decreasing steadily

**After Training:**
- Best model will be saved as `best_sep*.pt`
- Upload to S3
- Report results: separation, training time, final metrics

---

## Key Differences from Current Training

| Aspect | Current (1 GPU) | New (8 GPU) |
|--------|-----------------|-------------|
| Batch size per GPU | 1024 | 1024 |
| **Effective batch** | **1024** | **8192** ✅ |
| Training time | 3.44h (chr22 + full) | 5-8h (full only) |
| Expected accuracy | 96.6% | 99%+ |
| Separation | 12.40 | 13-15+ |

The key improvement is **effective batch size 8192** which provides more negative samples for InfoNCE contrastive learning!

---

## Files You'll Have

```
genocache-v4-8gpu/
├── GRCh38.fa              # Genome (3.1 GB)
├── GRCh38.fa.fai          # Index
├── SETUP_GUIDE.md         # Complete setup instructions
├── scripts/
│   ├── model.py           # Model architecture
│   ├── dataset.py         # Dataset with augmentation
│   ├── train_8gpu.py      # 8-GPU training script
│   └── validate.py        # Validation script
├── models/
│   └── fullgenome_best_sep12.4035_epoch30.pt  # Current best (reference)
├── validation/
│   └── test data
└── training.log           # Training progress
```

---

## Troubleshooting

**Problem:** CUDA out of memory  
**Solution:** Reduce batch size: `--batch-size 512` (effective batch = 4096)

**Problem:** Distributed training fails  
**Solution:** Check NCCL: `python3 -c "import torch; print(torch.distributed.is_nccl_available())"`

**Problem:** Slow data loading  
**Solution:** Increase workers: `--num-workers 16`

**Problem:** Loss not decreasing  
**Solution:** Check logs for errors, verify data loaded correctly

---

## Expected Timeline

```
00:00 - Setup environment (30 min)
00:30 - Download from S3 (20 min)
00:50 - Test single GPU (5 min)
01:00 - Launch 8-GPU training
01:00 - 09:00 - Training (5-8 hours)
09:00 - Upload to S3 (15 min)
09:15 - COMPLETE! 🎉
```

---

## What to Report Back

After training completes:

```
Training Complete! ✅

Results:
- Final separation: X.XX (vs 12.40 current)
- Final loss: X.XX
- Training time: X.X hours
- Best epoch: XX
- Effective batch size: 8192

Files uploaded to:
  s3://BUCKET/genocache-v4-8gpu/models/best_sep*.pt

Next steps:
- Validate on test data
- Compare with current model (96.6%)
- Expected improvement: 96.6% → 99%+
```

---

## Summary

**What you're doing:**
Training the same GenoCache V4 architecture with optimal InfoNCE batch size (8192 via 8 GPUs) to achieve 99%+ accuracy.

**Why it matters:**
- Current model: 96.6% (good, but not NeuralAligner level)
- Target model: 99%+ (matches state-of-the-art)
- Difference: Batch size 1024 → 8192

**Key commands:**
```bash
# Download
./download_from_s3.sh

# Setup
conda activate genocache && pip install dependencies

# Train
torchrun --nproc_per_node=8 scripts/train_8gpu.py --genome GRCh38.fa --batch-size 1024 --epochs 50

# Monitor
tail -f training.log

# Upload
aws s3 sync ./models/ s3://BUCKET/genocache-v4-8gpu/models/
```

**Expected outcome:**
GenoCache V4 with 99%+ accuracy, matching NeuralAligner performance! 🚀

---

**Good luck! You've got this! 💪**
