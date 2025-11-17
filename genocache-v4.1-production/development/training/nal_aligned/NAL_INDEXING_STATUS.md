# NAL-Aligned Indexing - IN PROGRESS

**Status:** Index building started  
**Date:** 2025-11-16 02:15 UTC  
**PID:** Check with `cat indexing.pid`

---

## Directory Structure

```
development/training/nal_aligned/
├── models/
│   └── genocache_nal.pt          ← Trained model (5.7MB, epoch 62)
├── indexes/
│   └── genocache_nal_stride32.index  ← Building... (ETA: ~45 min)
├── logs/
│   ├── training_nal_*.log        ← Training history
│   └── indexing_*.log            ← Indexing progress
├── encoder_nal.py                ← 128D encoder
├── augmentation_nal.py           ← Data augmentation
├── train_nal.py                  ← Training script
├── build_index_nal.py            ← Index builder (following NAL)
├── START_INDEXING.sh             ← Indexing launcher
└── README.md                     ← Workflow docs
```

**Everything organized in one place!** ✅

---

## Indexing Configuration (Following NAL Paper)

Based on Section 3.2 and A.3 of NeuralAligner paper:

### Index Parameters

```python
# Following NAL exactly:
index_type = "IVFPQ"          # Inverted File + Product Quantization
seed_length = 512             # NAL uses 256-512bp
stride = 32                   # NAL: ≤ seed_len/8 (≤64 for 512bp)
dimension = 128               # NAL uses 128D
compression = "PQ16×8"        # 16 subquantizers × 8 bits
nlist = sqrt(num_vectors)     # Auto-calculated
nprobe = 8                    # NAL recommends 8-32
distance = "inner_product"    # NAL uses inner product
use_projection = False        # CRITICAL: Use encoder output directly!
```

### Why These Parameters?

**Seed length (512bp):**
- Longer seeds = more specific, fewer false matches
- NAL shows 512bp gives best recall@1 (95.57% vs 94.16% for 256bp)
- Table 4 in paper demonstrates this clearly

**Stride (32bp):**
- NAL: "stride should be no larger than 1/8 of seed length"
- For 512bp: max stride = 64bp
- We use 32bp for better coverage (conservative)
- Translation continuity allows sparse indexing

**IVFPQ compression:**
- Reduces 128D×4bytes (512 bytes) → 16 bytes
- GRCh38 with stride=32: ~90M vectors × 16 bytes = ~1.5GB
- Fits in GPU memory for fast search
- NAL paper: "only 2.1 GB for GRCh38" with PQ16×8

**nprobe = 8:**
- Lower = faster, slightly less recall
- NAL: "between 8 and 32"
- We use 8 for speed (can increase to 32 for accuracy)

**NO projection head:**
- Projection is ONLY for training (prevents collapse)
- During indexing/inference: Use encoder output directly
- This is critical! NAL paper Section 3.1 explicitly states this

---

## Expected Output

### Estimated Statistics

```
Genome: GRCh38.p14 (2.65 Gbp)
Chromosomes: 676 sequences

Seeds extracted:
  • Total: ~82.8M seeds
  • Calculation: 2,649,480,772 bp / 32 stride = 82,796,274
  • After filtering (high-N): ~80M-85M

Index size:
  • Uncompressed: 82M × 128D × 4 bytes = 42 GB
  • With PQ16×8: 82M × 16 bytes = 1.3 GB ✅
  
Search time (per seed):
  • With nprobe=8: ~2-5 μs (microseconds)
  • NAL paper Table 7: 2.5 μs per seed query
```

### Build Time

NAL paper (Section A.3):
> "Indexing the entire human genome on a single NVIDIA H20 GPU takes  
> between 20 and 50 minutes, depending on the configuration."

**Our estimate:** ~45 minutes on H100 GPU
- Loading genome: ~2 min
- Extracting seeds: ~10 min
- Encoding seeds: ~25 min (80M seeds, batch=2048)
- Building FAISS index: ~8 min
- **Total:** ~45 minutes

---

## Monitoring

### Quick Check
```bash
/home/nebius/genocache/CHECK_INDEXING.sh
```

### Live Monitoring
```bash
watch -n 10 /home/nebius/genocache/CHECK_INDEXING.sh
```

### Check Process
```bash
ps aux | grep build_index_nal
```

### Check Logs
```bash
tail -f development/training/nal_aligned/indexing_run.log
```

---

## What Happens During Indexing

### Phase 1: Load Model & Genome (2 min)
```
✓ Load trained NAL model (epoch 62, val_loss=0.19)
✓ Load GRCh38.fa (2.65 Gbp, 676 sequences)
✓ Filter out high-N regions (>10% N)
```

### Phase 2: Extract Seeds (10 min)
```
For each chromosome:
  Slide 512bp window with stride=32
  Skip if >10% N
  Store: (chr_name, position, seed_sequence)
  
Output: ~82M seeds + positions
```

### Phase 3: Encode Seeds (25 min)
```
For each batch of 2048 seeds:
  DNA → token IDs (A=0, C=1, G=2, T=3, N=4)
  Token embedding → 4× conv → pool → L2 norm
  Output: 128D embedding (NO projection!)
  
Output: 82M × 128D = 10.5 GB of embeddings (FP32)
```

### Phase 4: Build FAISS Index (8 min)
```
1. Calculate nlist = sqrt(82M) ≈ 9,055 clusters
2. Train quantizer on 256 × nlist samples
3. Normalize embeddings (L2 norm)
4. Add all vectors with PQ16×8 compression
5. Set nprobe = 8
6. Save index + position mapping
```

---

## After Indexing Completes

### Files Created

```
indexes/genocache_nal_stride32.index          (~1.3 GB)
indexes/genocache_nal_stride32.positions.npz  (~1.5 GB)
```

### Next Step: Test on GIAB Data

```bash
# Update alignment pipeline to use new index
python3 genocache_align_minimap2.py \
    --reads validation/data/giab_hg002_100reads.fastq \
    --output validation/results/giab_nal.sam \
    --reference /home/nebius/genocache/GRCh38.fa \
    --model development/training/nal_aligned/models/genocache_nal.pt \
    --index development/training/nal_aligned/indexes/genocache_nal_stride32.index
```

**Expected result:** 85-90% mapping rate (vs 32% current, 94% minimap2)

---

## NAL Paper References

### Section 3.2: Embedding-Based Genome Indexing

> "We further take advantage of the translation continuity of the embeddings  
> to build a sparse index for faster search. Instead of indexing every position  
> in the genome with a sliding window of step size 1, we use a window of size k.  
> Each seed embedding then represents k consecutive bases, and for each window,  
> we index only the embedding of the seed centered at the middle base."

### Section A.3: Indexing Configuration

> "We build the index using Faiss (Johnson et al., 2019), with inverted file  
> and product quantization (IVFPQ)... With product quantization and an indexing  
> stride of 32, the GRCh38.p14 primary assembly requires only 2.1 GB of GPU  
> memory, which is feasible on most modern GPUs."

**Guidelines:**
- nlist = sqrt(num_vectors)
- stride ≤ seed_len/8
- nprobe = 8-32
- K = 32 (top-K neighbors)
- Distance: inner product

---

## Troubleshooting

### If indexing fails:

1. **Check logs:**
   ```bash
   tail -100 development/training/nal_aligned/indexing_run.log
   ```

2. **Check GPU memory:**
   ```bash
   nvidia-smi
   ```
   
3. **Restart with smaller batch:**
   ```bash
   # Edit START_INDEXING.sh: change --batch-size 2048 to 1024
   ./START_INDEXING.sh
   ```

### If index is too large:

- Increase stride to 64 (reduces vectors by 50%)
- Use PQ32×4 instead of PQ16×8 (smaller, slightly less accurate)

---

## Summary

**Status:** ✅ Index building in progress  
**Configuration:** Exact NAL protocol (Section A.3)  
**Output:** IVFPQ index with PQ16×8 compression  
**Size:** ~1.3 GB (compressed from 42 GB)  
**Coverage:** ~82M seeds across 2.65 Gbp genome  
**ETA:** ~45 minutes  
**Next:** Test on GIAB, expect 85-90% mapping rate
