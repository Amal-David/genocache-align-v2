# NAL-Aligned Training Pipeline

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                    TRAINING LOOP (Per Batch)                     │
└─────────────────────────────────────────────────────────────────┘

1. DATA GENERATION (CPU - Parallel):
   ┌──────────────────────────────────────────────┐
   │ genome → extract 512bp random position       │
   │       ↓                                       │
   │ anchor (original)  →  positive (augmented)   │
   │                                               │
   │ Augmentation:                                 │
   │  • Add errors: U[0.01, 0.1] rate             │
   │  • Shift: ±51bp                              │
   │  • RC flip: 50% each independently           │
   │                                               │
   │ Output: (anchors[1024], positives[1024])     │
   └──────────────────────────────────────────────┘

2. ENCODING (GPU):
   ┌──────────────────────────────────────────────┐
   │ DNA string → Token IDs (A=0,C=1,G=2,T=3,N=4) │
   │           ↓                                   │
   │ Token Embedding (5 → 128D)                   │
   │           ↓                                   │
   │ 4× Bidirectional Conv (kernel=7)             │
   │           ↓                                   │
   │ Average Pool → L2 Normalize → 128D           │
   │           ↓                                   │
   │ Projection Head (training only) → 128D       │
   │                                               │
   │ Output: embeddings[1024, 128]                │
   └──────────────────────────────────────────────┘

3. LOSS COMPUTATION (GPU):
   ┌──────────────────────────────────────────────┐
   │ InfoNCE Loss:                                 │
   │                                               │
   │ similarity = anchor · positive^T / temp       │
   │ loss = -mean(log(exp(pos) / sum(exp(all))))  │
   │                                               │
   │ Effect:                                       │
   │  • Pull anchor-positive pairs together        │
   │  • Push all negatives (other samples) apart   │
   └──────────────────────────────────────────────┘

4. OPTIMIZATION (GPU):
   ┌──────────────────────────────────────────────┐
   │ AdamW optimizer updates weights               │
   │ lr=1e-3, weight_decay=0.1                    │
   └──────────────────────────────────────────────┘

128 batches = 1 epoch
~100 epochs = 13M sample pairs (NAL spec)
```

## Files

- `encoder_nal.py` - 128D encoder (493K params)
- `augmentation_nal.py` - Data augmentation
- `train_nal.py` - Training loop

## Usage

```bash
python3 train_nal.py \
    --genome /home/nebius/genocache/GRCh38.fa \
    --output ../../models/genocache_nal.pt \
    --device cuda
```
