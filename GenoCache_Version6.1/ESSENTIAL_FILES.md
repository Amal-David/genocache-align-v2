# Essential Files for GenoCache v6.1

## Location: `/home/nebius/genocache/GenoCache_Version6.1/`

## Core Training & Indexing Files

### 1. Training
```
/home/nebius/genocache/GenoCache_Version6.1/training/train_full_genome.py
```
- **Purpose:** Full genome training (16,000 batches)
- **Status:** ✅ Ready to use
- **Features:** GPU auto-detection, multi-GPU support, curriculum learning

### 2. Encoder Architecture
```
/home/nebius/genocache/GenoCache_Version6.1/training/encoder_nal.py
```
- **Purpose:** Neural encoder model (128D embeddings)
- **Status:** ✅ Proven (76.5% accuracy on Chr22)
- **Size:** 34 lines (simplified, core only)

### 3. Index Building (NEED TO CREATE)
```
/home/nebius/genocache/GenoCache_Version6.1/indexing/build_full_genome_index.py
```
- **Purpose:** Build FAISS index from trained model
- **Status:** ⏳ To be created based on Chr22 success

### 4. Alignment Pipeline (NEED TO COPY)
```
/home/nebius/genocache/GenoCache_Version6.1/alignment/seeding.py
/home/nebius/genocache/GenoCache_Version6.1/alignment/chaining.py
/home/nebius/genocache/GenoCache_Version6.1/alignment/alignment.py
```
- **Purpose:** Complete alignment pipeline
- **Status:** ⏳ To be copied from successful Chr22 version

## Reference Files (Where It Worked)

### Chr22 Success Files (76.5% accuracy)
```
/home/nebius/genocache/genocache-v4.1-production/development/training/nal_chr22_PRODUCTION_BACKUP/
```

Contains:
- `training/train_chr22_focused.py` - Training that got 76.5%
- `training/encoder_nal.py` - Same encoder (34 lines)
- `indexing/build_chr22_index.py` - Index building
- `alignment/seeding_nal.py` - Seeding phase
- `alignment/chaining_nal.py` - Chaining phase
- `documentation/` - All documentation

### Working Chr22 Models & Indexes
```
/home/nebius/genocache/genocache-v4.1-production/development/training/nal_single_chr/
```

Contains:
- `models/` - Trained models (5.5 MB each)
- `indexes/` - FAISS indexes (1.1 GB)
- `results/` - Test results showing 76.5% accuracy

## Files You DON'T Need

❌ `/home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned/`
- This was the ORIGINAL attempt (43% accuracy, 88.89% false positives)
- DO NOT USE these files
- Chr22 version is the one that works!

## Quick Copy Commands

### To backup the working Chr22 version:
```bash
# Create backup
tar -czf geocache_chr22_success.tar.gz \
  /home/nebius/genocache/genocache-v4.1-production/development/training/nal_chr22_PRODUCTION_BACKUP/

# Check size
ls -lh geocache_chr22_success.tar.gz
```

### To backup v6.1 (after training):
```bash
tar -czf geocache_v6.1_full_genome.tar.gz \
  /home/nebius/genocache/GenoCache_Version6.1/

ls -lh geocache_v6.1_full_genome.tar.gz
```

## What's Next

1. ✅ Training files ready (train_full_genome.py)
2. ⏳ Create index building script
3. ⏳ Copy alignment pipeline from Chr22
4. ⏳ Test full pipeline
5. ⏳ Benchmark against baselines

## File Sizes

- Training script: ~15 KB
- Encoder: ~2 KB
- Trained model: ~5-10 MB
- Full genome index: ~20-30 GB
- Complete v6.1: ~30-40 GB (after training)

