# Development Scripts - Training & Indexing

Scripts used to train the model and build the production index.

---

## Models in Production

**Current Model (in use):**
- File: `../models/genocache_model.pt`
- Source: `fullgenome_best_sep12.3547_epoch8.pt`
- Separation Score: 12.3547
- Status: ✅ Validated (37%→87% with EXTEND)

**Latest Model (newer, better):**
- File: `../models/genocache_model_latest.pt`
- Source: `fullgenome_best_sep12.4035_epoch30.pt`
- Separation Score: 12.4035 (better!)
- Date: Nov 13, 2024
- Status: ⚠️ Not yet validated with EXTEND

**Recommendation:** Test `genocache_model_latest.pt` to see if it improves accuracy further!

---

## Training Scripts

### Main Training Script (Used for Production Model)

**`training/train_fullgenome.py`** ⭐
- **What it does:** Trains model on full GRCh38 genome (all 25 chromosomes)
- **Output:** `fullgenome_best_*.pt` checkpoints
- **Training method:** Contrastive learning (InfoNCE loss)
- **Data augmentation:** See `augmentation.py`

**Usage:**
```bash
cd training
python train_fullgenome.py
```

**Key parameters:**
- Batch size: 1024
- Learning rate: 1e-3 with scheduler
- Epochs: 30
- Embedding dim: 128D
- Seed length: 512bp

### Supporting Training Scripts

**`training/augmentation.py`**
- Data augmentation pipeline
- Creates positive pairs with noise/shifts
- Error rates: 1-10% (mimics sequencing errors)
- Shift range: ±51bp (10% of seed length)

**`training/train_chr22.py`**
- Proof-of-concept training on chr22 only
- Faster iteration for testing
- Used to validate architecture before full genome

**`training/train_curriculum.py`**
- Curriculum learning approach
- Starts with easy examples, gradually harder
- Experimental (not used for production)

**`training/train_simple.py`**
- Minimal training script
- Good for understanding basics
- Not used for production

**`training/generate_training_data.py`**
- Generates pre-computed training datasets
- Alternative to on-the-fly generation
- Not used for current model (we use on-the-fly)

---

## Indexing Scripts

### Production Index Builder

**`indexing/build_production_index.py`** ⭐
- **What it does:** Builds FAISS IVFPQ index from model
- **Output:** `genocache_v4_production.index` + metadata
- **Index type:** IVFPQ (Inverted File + Product Quantization)
- **Compression:** PQ16x8 (16 subvectors, 8-bit codes)

**Usage:**
```bash
cd indexing
python build_production_index.py \
    --model ../models/genocache_model.pt \
    --reference /path/to/GRCh38.fa \
    --output ../indexes/genocache_v4_production
```

**What it does:**
1. Loads trained model
2. Extracts 512bp windows from genome (stride=32)
3. Encodes each window → 128D embedding
4. Builds FAISS index with ~91.8M vectors
5. Saves index + metadata (positions, chr names)

**Time:** ~2-4 hours on GPU

### Alternative Index Builder

**`indexing/build_indexes.py`**
- General-purpose index builder
- Supports multiple index types (Flat, IVF, HNSW)
- Used for experiments
- Not used for production (use `build_production_index.py`)

---

## Training Pipeline (Full Workflow)

### Step 1: Prepare Data

**Download GRCh38:**
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa
```

### Step 2: Train Model

```bash
cd development/training

# Train on full genome (recommended)
python train_fullgenome.py

# Or test on chr22 first (faster)
python train_chr22.py
```

**Output:**
- Checkpoints saved to `../../models/checkpoints/`
- Best model selected by separation score (pos_sim - neg_sim)
- Training takes ~12-24 hours on single GPU

### Step 3: Build Index

```bash
cd development/indexing

# Build production index
python build_production_index.py \
    --model ../../models/checkpoints/fullgenome_best_*.pt \
    --reference /path/to/GRCh38.fa \
    --output ../../indexes/genocache_v4_production
```

**Output:**
- `genocache_v4_production.index` (2.1GB)
- `genocache_v4_production.metadata.pkl` (4.8GB)
- Takes ~2-4 hours on GPU

### Step 4: Validate

```bash
cd ../..

# Test alignment with new model
python genocache_align.py \
    --reads test_reads.fastq \
    --output test.sam \
    --reference /path/to/GRCh38.fa \
    --model models/checkpoints/fullgenome_best_*.pt \
    --index indexes/genocache_v4_production.index \
    --metadata indexes/genocache_v4_production.metadata.pkl

# Run validation
python scripts/validate_extend_fix.py
```

---

## Training Details

### Architecture (models/encoder.py)

```python
GenoCacheEncoder(
    emb_dim=128,              # Output embedding size
    seed_len=512,             # Input sequence length
    vocab_size=5,             # A, C, G, T, N
    hidden_dims=[64, 128, 256],  # Multi-scale features
    num_attention_layers=2,   # Lightweight attention
    dropout=0.1
)
```

**Parameters:** 1.2M  
**Architecture:** Multi-scale CNN + Lightweight Attention

### Data Augmentation

**From `augmentation.py`:**

```python
# For each training sample:
1. Extract 512bp window from genome
2. Create positive pair:
   - Add sequencing errors (1-10% error rate)
     - Substitutions, insertions, deletions
   - Shift position (±51bp, 10% of seed length)
   - Pad/trim to 512bp
3. Encode both → 128D embeddings
4. InfoNCE loss: pull positives together, push negatives apart
```

**Why this works:**
- Robustness to sequencing errors (trained with 1-10% errors)
- Translation continuity (trained with position shifts)
- Longer seeds (512bp) more informative than traditional 15-19bp

### Loss Function

**InfoNCE (Contrastive Learning):**

```python
# For batch of N samples:
similarity = anchor @ positive.T / temperature

# Positive similarity (diagonal)
pos_sim = diagonal(similarity)

# Loss: maximize pos_sim, minimize neg_sim
loss = -mean(pos_sim - logsumexp(similarity, dim=1))
```

**Separation score:** `pos_sim - neg_sim`
- Higher = better model
- Production model: 12.3547 (epoch 8)
- Latest model: 12.4035 (epoch 30)

### Training Hyperparameters

```python
# From train_fullgenome.py:
batch_size = 1024
learning_rate = 1e-3
optimizer = AdamW(lr=1e-3, weight_decay=0.1)
scheduler = ReduceLROnPlateau(patience=4, factor=0.2)
temperature = 0.07  # InfoNCE temperature
epochs = 30
batches_per_epoch = 500
```

**Training time:** ~30 epochs × 500 batches × 2 sec/batch = ~8-12 hours

---

## Indexing Details

### Index Type: FAISS IVFPQ

**Parameters:**
```python
# From build_production_index.py:
nlist = int(sqrt(n_vectors))  # ~9600 clusters
nprobe = 16                    # Search 16 clusters
m = 16                         # 16 subvectors
nbits = 8                      # 8-bit codes
```

**Compression:**
- Original: 128D × 4 bytes = 512 bytes/vector
- Compressed: 16 bytes/vector (PQ16x8)
- Ratio: 32:1 compression

**Index size:**
- Vectors: 91.8M
- Size: 91.8M × 16 bytes = ~1.5GB (actual: 2.1GB with overhead)
- Metadata: 91.8M × 52 bytes = ~4.8GB (positions + chr names)

### Indexing Strategy

**Sliding window:**
- Window size: 512bp
- Stride: 32bp (1/16 of window)
- Why: Balance coverage vs index size
  - Stride=1: Perfect coverage, 32× more vectors
  - Stride=32: Good coverage (translation continuity helps), reasonable size

**Coverage:**
- GRCh38: 3.1 billion bp
- Windows: 3.1B / 32 = ~97M positions
- Actual: 91.8M (some regions skipped, e.g., N-rich)

---

## Reproducing Production Model

**Complete workflow:**

```bash
# 1. Download reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/.../GRCh38.fa.gz
gunzip GRCh38.fa.gz

# 2. Train model
cd development/training
python train_fullgenome.py
# → Produces: fullgenome_best_*.pt

# 3. Build index
cd ../indexing
python build_production_index.py \
    --model ../../models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt \
    --reference /path/to/GRCh38.fa \
    --output ../../indexes/genocache_v4_production_new

# 4. Test
cd ../..
python genocache_align.py \
    --model models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt \
    --index indexes/genocache_v4_production_new.index \
    --metadata indexes/genocache_v4_production_new.metadata.pkl \
    --reads test.fastq \
    --output test.sam \
    --reference /path/to/GRCh38.fa

# 5. Validate
python scripts/validate_extend_fix.py
```

---

## Files Included

**Training:**
- ✅ `training/train_fullgenome.py` - Main training script ⭐
- ✅ `training/augmentation.py` - Data augmentation
- ✅ `training/train_chr22.py` - Chr22-only training
- ✅ `training/train_curriculum.py` - Curriculum learning
- ✅ `training/train_simple.py` - Minimal example
- ✅ `training/generate_training_data.py` - Data generation
- ✅ `training/train_10M.py` - 10M sample training
- ✅ `training/generate_training_data_10M.py` - 10M data gen

**Indexing:**
- ✅ `indexing/build_production_index.py` - Production index builder ⭐
- ✅ `indexing/build_indexes.py` - General index builder

---

## Next Steps

### Test Latest Model

The latest model (`fullgenome_best_sep12.4035_epoch30.pt`) has better separation score (12.4035 vs 12.3547).

**To test:**
1. Copy to production: `cp ../models/genocache_model_latest.pt ../models/genocache_model.pt`
2. Run validation: `python ../scripts/validate_extend_fix.py`
3. Check if accuracy improves beyond 87.5%

### Retrain from Scratch (If Needed)

If you need to retrain:
```bash
cd development/training
python train_fullgenome.py
```

Then rebuild index:
```bash
cd ../indexing
python build_production_index.py --model <new_model> --reference GRCh38.fa
```

### Improve Training (Ideas)

1. **Longer training:** Try 50+ epochs
2. **Larger batch size:** Try 2048 (if GPU memory allows)
3. **Hard negative mining:** Sample harder negatives
4. **Multi-GPU:** Distribute training across GPUs
5. **Different augmentation:** Try different error rates/shifts

---

## Notes

### Model Selection

**Separation score (pos_sim - neg_sim):**
- Higher = better discrimination
- Production: 12.3547
- Latest: 12.4035 ← Test this!

**When to stop training:**
- Separation score plateaus
- Validation performance stops improving
- Overfitting detected

### Index Updates

**When to rebuild index:**
- New model trained
- Reference genome updated
- Different index parameters needed

**Index rebuild takes ~2-4 hours** - plan accordingly!

---

**Last Updated:** 2025-11-15  
**Version:** 4.1.0  
**Status:** Production Development Scripts
