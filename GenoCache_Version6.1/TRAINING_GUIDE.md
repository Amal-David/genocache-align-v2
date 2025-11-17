# GenoCache v6.1 - Training Guide

## Quick Start

### Start Training Now

```bash
cd /home/nebius/genocache/GenoCache_Version6.1/training
python3 train_full_genome.py
```

This will train for **24,000 batches** (~4-6 hours).

## Training Parameters

### What's Different from Chr22 (v4.1)

| Parameter | Chr22 (v4.1) | Full Genome (v6.1) | Why Changed |
|-----------|--------------|-------------------|-------------|
| Chromosomes | 1 (Chr22) | 24 (All) | Scale to full genome |
| Batches | 8,000 | **24,000** | 3x more for larger genome |
| Curriculum | 5→15% (4 stages) | 5→15% (4 stages) | Same proven approach |
| Training time | 1.5 hours | 4-6 hours | 3x genome size |
| Expected accuracy | 76.5% | 75-85% | Slightly lower (harder) |

### Curriculum Learning Schedule

```
Batches 0-3,000:       5% error rate (EASY - learn clean patterns)
Batches 3,000-8,000:   8% error rate
Batches 8,000-16,000:  12% error rate
Batches 16,000-24,000: 15% error rate (HARD - handle ONT errors)
```

### What Stays the Same (Proven in Chr22)

✅ **Seed length:** 512bp  
✅ **Architecture:** CNN-based, 4 layers, 128D embeddings  
✅ **Error injection:** Sub/Ins/Del independently  
✅ **Loss function:** InfoNCE contrastive  
✅ **Batch size:** 32  
✅ **Learning rate:** 1e-4

## Monitoring Training

Training will print progress every 100 batches:

```
Batch 100/24,000 | Error: 5.0% | Loss: 0.8234 | Acc: 72.3% | ETA: 4.2h
Batch 200/24,000 | Error: 5.0% | Loss: 0.6541 | Acc: 78.9% | ETA: 4.1h
...
```

Checkpoints saved every 2,000 batches:
- `models/genocache_nal_batch2000.pt`
- `models/genocache_nal_batch4000.pt`
- ...
- `models/genocache_nal_fullgenome_final.pt`

## Estimated Timeline

```
Hour 0:     Start training (5% errors)
Hour 0.5:   Batch 3,000 - switch to 8% errors
Hour 1.3:   Batch 8,000 - switch to 12% errors
Hour 2.7:   Batch 16,000 - switch to 15% errors
Hour 4-6:   Complete! (24,000 batches done)
```

## After Training

### 1. Build Index (~2-3 hours)
```bash
cd ../indexing
python3 build_full_genome_index.py
```

### 2. Test Alignment
```bash
cd ../alignment
python3 test_alignment.py
```

### 3. Expected Results

Based on Chr22 success (76.5%), we expect:
- **Accuracy:** 75-85% on synthetic reads
- **vs minimap2:** Competitive or better
- **Speed:** 5-10x faster than minimap2

## Resource Requirements

- **GPU:** NVIDIA H100 (or similar)
- **Memory:** ~40 GB RAM
- **Disk:** ~50 GB free space
- **Time:** 4-6 hours training + 2-3 hours indexing

## Troubleshooting

### Training Too Slow?
- Check GPU utilization: `nvidia-smi`
- Reduce batch_size from 32 to 16
- Or train on fewer chromosomes first

### Out of Memory?
- Reduce batch_size: Edit `train_full_genome.py` line with `batch_size=32`
- Change to `batch_size=16`

### Want to Resume?
- Training saves checkpoints every 2000 batches
- Modify script to load from checkpoint
