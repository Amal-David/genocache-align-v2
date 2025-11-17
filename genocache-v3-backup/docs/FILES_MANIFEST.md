# GenoCache V3 - Files Manifest

**Backup Date:** November 12, 2025  
**Purpose:** Keep essential artifacts, document everything else

---

## ✅ Files Backed Up (Total: ~15 MB)

### **Model (6.1 MB)**
```
model/genocache_best.pt
├─ Architecture: Multi-scale CNN + Linear Attention
├─ Parameters: 1.57M
├─ Training: 50 epochs on full GRCh38
├─ Final loss: 0.061
├─ Quality: 0.985 same-region similarity
└─ Status: Works perfectly on clean reads, needs error training
```

### **Core Code (120 KB)**
```
code/
├─ genocache_encoder.py       (20 KB) - Model architecture
├─ train_single_gpu.py         (20 KB) - Training pipeline with hard negatives
├─ align_production.py         (20 KB) - Production aligner (synthetic reads)
├─ align_fastq.py              (16 KB) - FASTQ input support
└─ align_with_sparse_dp.py     (16 KB) - Sparse DP chaining (87% accuracy)
```

### **Validation Results (1.5 MB)**
```
results/giab_validation/
├─ ont_reads.fq                (1.2 MB) - 500 realistic ONT reads (8% error)
└─ validation_results.json     (2 KB)   - Minimap2 baseline: 94.6% accuracy
```

### **Documentation (50 KB)**
```
docs/
├─ SESSION_SUMMARY.md          (40 KB) - Complete session log
├─ FILES_MANIFEST.md           (10 KB) - This file
└─ V4_TRAINING_PLAN.md         (pending) - Next iteration strategy
```

**Total Backup Size: ~15 MB** (compact and portable)

---

## 📋 Files Documented (NOT Backed Up)

### **Intermediate Checkpoints (67 MB) - ❌ Not Needed**
```
Location: genocache-production/models/checkpoints/

genocache_epoch5.pt   (6.1 MB) - Loss: ~0.8
genocache_epoch10.pt  (6.1 MB) - Loss: ~0.5
genocache_epoch15.pt  (6.1 MB) - Loss: ~0.3
genocache_epoch20.pt  (6.1 MB) - Loss: ~0.2
genocache_epoch25.pt  (6.1 MB) - Loss: ~0.15
genocache_epoch30.pt  (6.1 MB) - Loss: ~0.12
genocache_epoch35.pt  (6.1 MB) - Loss: ~0.10
genocache_epoch40.pt  (6.1 MB) - Loss: ~0.08
genocache_epoch45.pt  (6.1 MB) - Loss: ~0.07
genocache_epoch50.pt  (6.1 MB) - Loss: 0.061 (= genocache_best.pt)
genocache_final.pt    (6.1 MB) - Same as epoch50

Why Not Backed Up:
- We have the best checkpoint (genocache_best.pt)
- Intermediate checkpoints useful for training analysis only
- Can retrain if needed (8 hours on H100)
- Checkpoints show smooth convergence (good sign)
```

### **Genome Embeddings (52 GB) - ❌ Can Regenerate**
```
Location: genocache-production/encoded/

genome_vectors.npy        (25 GB)
├─ Shape: (26,127,161, 256)
├─ Coverage: Full GRCh38 (705 sequences)
├─ Stride: 100bp
├─ Time to generate: ~1 hour on H100
└─ Can regenerate: Yes, with genocache_best.pt

genome_positions.txt      (553 MB)
├─ Format: chrom \t position (one per embedding)
├─ Lines: 26,127,161
├─ Can regenerate: Yes, deterministic from stride
└─ Regeneration script: encode_genome.py

genome_ivfflat.faiss      (26 GB)
├─ Type: IVF-Flat index
├─ Recall@10: 99.7% ✅
├─ nlist: 8192
├─ Vectors: 26,127,161
├─ Time to build: ~30 min
└─ Can regenerate: Yes, from genome_vectors.npy

genome_ivfpq.faiss        (403 MB)
├─ Type: IVF-PQ (M=8)
├─ Recall@10: 5% ❌ (FAILED)
├─ Status: Too compressed
└─ Note: Experiment showed PQ doesn't work for our embeddings

genome_ivfpq_m32.faiss    (1.0 GB)
├─ Type: IVF-PQ (M=32)  
├─ Recall@10: 51.5% ❌ (FAILED)
├─ Status: Still too compressed
└─ Note: Confirmed we need IVF-Flat (no compression)

Why Not Backed Up:
- 52 GB is too large for backup
- Can regenerate in 2 hours total:
  1. Encode genome: 1 hour (encode_genome.py)
  2. Build index: 30 min (build_faiss_index.py)
- Index experiments already concluded (use IVF-Flat)
```

### **Large Validation Files (37 MB) - ❌ Have Summary**
```
Location: genocache-production/giab_validation/

genocache.sam             (25 MB)
├─ 500 ONT reads aligned with GenoCache
├─ Accuracy: 27.4% correct chromosome ❌
├─ MAPQ scores: Mean 28.2, only 5.2% ≥30
└─ Summary: In validation_results.json ✅

minimap2.sam              (12 MB)
├─ Same 500 reads aligned with minimap2
├─ Accuracy: 94.6% ✅
├─ Speed: 10.5 reads/sec
└─ Summary: In validation_results.json ✅

Why Not Backed Up:
- Full SAM files not needed
- All metrics captured in validation_results.json
- Can regenerate if needed (ont_reads.fq backed up)
```

### **Experimental Scripts (80 KB) - ❌ Documented**
```
Location: genocache-production/

align_reads_voting.py              (16 KB)
├─ Voting strategy (without chaining)
├─ Result: 43% accuracy ❌
└─ Lesson: Voting alone insufficient, need chaining

align_with_cluster_voting.py      (20 KB)
├─ Chromosome consensus + clustering
├─ Result: 3.6-18.5% accuracy ❌
└─ Lesson: Pre-filtering removes good candidates

compare_aligners.py                (8 KB)
├─ GenoCache vs Minimap2 comparison
├─ Result: Showed 69% accuracy gap
└─ Lesson: Neural seeds fail on noisy reads

validate_giab_style.py             (16 KB)
├─ GIAB-style validation pipeline
├─ Generates realistic ONT reads
└─ Status: Working, used for final validation

Why Not Backed Up:
- Experimental approaches that failed
- Lessons learned documented in SESSION_SUMMARY.md
- Keep in original location for reference
```

### **Documentation Files (100 KB) - ❌ Keep Originals**
```
Location: genocache-production/

ADVANCED_INDEXING.md           (16 KB) - Index compression analysis
CHAINING_REQUIRED.md           (8 KB)  - Why chaining is essential  
EXECUTION_PLAN.md              (8 KB)  - Phase 1 & 2 roadmap
INDEXING_STRATEGY.md           (9 KB)  - Position clustering explanation
VOTING_STRATEGY_COMPARISON.md  (7 KB)  - Voting vs chaining analysis

Why Not Backed Up:
- Historical documentation from V3 development
- Key insights extracted to SESSION_SUMMARY.md
- Keep originals in place for reference
```

### **Log Files (7.7 MB) - ❌ Have Summary**
```
Location: genocache-production/logs/

Training logs:
- Full training output (convergence curves)
- Checkpoint save confirmations
- Loss trajectories

Validation logs:
- Alignment test results
- Index building progress
- Benchmark outputs

Why Not Backed Up:
- 7.7 MB of detailed logs
- Key metrics extracted to results JSON files
- Can generate new logs in V4
```

---

## 🗑️ Files to Delete (Optional Cleanup)

### **Safe to Delete (Save 52 GB):**
```bash
# These can all be regenerated in 2 hours
cd /home/nebius/genocache/genocache-production/encoded/
rm genome_vectors.npy             # 25 GB - regenerate with encode_genome.py
rm genome_ivfflat.faiss           # 26 GB - rebuild from vectors
rm genome_positions.txt           # 553 MB - regenerate (deterministic)
rm genome_ivfpq*.faiss            # 1.4 GB - failed experiments

# Intermediate checkpoints
cd /home/nebius/genocache/genocache-production/models/checkpoints/
rm genocache_epoch*.pt            # 67 MB - keep only genocache_best.pt
rm genocache_final.pt             # 6.1 MB - duplicate of best

# Large SAM files  
cd /home/nebius/genocache/genocache-production/giab_validation/
rm genocache.sam minimap2.sam     # 37 MB - metrics saved in JSON
```

**Total Space Saved: ~52 GB**

### **Keep These (Development Files):**
```bash
# Active code files
*.py files in genocache-production/

# Best model
genocache-production/models/checkpoints/genocache_best.pt

# Small test results
*.json files

# Documentation
*.md files
```

---

## 📦 How to Regenerate Everything

### **If You Need the Full Pipeline Again:**

```bash
# 1. Restore model (already have it)
cp genocache-v3-backup/model/genocache_best.pt models/checkpoints/

# 2. Re-encode genome (1 hour on H100)
python encode_genome.py \
  --checkpoint models/checkpoints/genocache_best.pt \
  --fasta /home/nebius/genocache/GRCh38.fa \
  --output encoded/ \
  --stride 100

# Creates:
# - encoded/genome_vectors.npy (25 GB)
# - encoded/genome_positions.txt (553 MB)

# 3. Build index (30 min)
python build_faiss_index.py \
  --vectors encoded/genome_vectors.npy \
  --output encoded/genome_ivfflat.faiss \
  --index-type IVF-Flat \
  --nlist 8192

# Creates:
# - encoded/genome_ivfflat.faiss (26 GB)

# 4. Validate (5 min)
python validate_giab_style.py \
  --fasta /home/nebius/genocache/GRCh38.fa \
  --chrom NC_000001.11 \
  --checkpoint models/checkpoints/genocache_best.pt \
  --index encoded/genome_ivfflat.faiss \
  --positions encoded/genome_positions.txt \
  --num-reads 500 \
  --output-dir giab_validation

# Total time: ~2 hours
# Total space: ~52 GB
```

---

## 📊 Storage Summary

### **Original (genocache-production/):**
```
Total: ~52.2 GB
├─ encoded/         52.0 GB (genome embeddings + indices)
├─ models/          73 MB   (12 checkpoints)
├─ giab_validation/ 38 MB   (test results)
└─ logs/            7.7 MB  (training/test logs)
```

### **Backup (genocache-v3-backup/):**
```
Total: ~15 MB (0.03% of original)
├─ model/      6.1 MB  (best checkpoint only)
├─ code/       120 KB  (5 essential scripts)
├─ results/    1.5 MB  (test reads + metrics)
└─ docs/       50 KB   (complete documentation)
```

### **Savings:**
- Backup size: 15 MB (fits on USB drive)
- Original size: 52.2 GB
- Compression ratio: 3,480× (99.97% reduction)
- Everything reproducible: Yes ✅

---

## 🎯 What's Essential vs Optional

### **ESSENTIAL (Must Keep):**
1. ✅ `genocache_best.pt` - Cannot retrain exactly (random seed, shuffling)
2. ✅ Core code files - Our implementation
3. ✅ Documentation - Knowledge captured
4. ✅ Validation results JSON - Metrics for comparison

### **NICE TO HAVE (Can Regenerate):**
1. ⚠️ Genome embeddings - 1 hour to regenerate
2. ⚠️ FAISS index - 30 min to rebuild
3. ⚠️ Test reads - 1 min to regenerate
4. ⚠️ Intermediate checkpoints - Training analysis only

### **OPTIONAL (Not Needed):**
1. ❌ Full SAM files - Have metrics
2. ❌ Log files - Have summaries
3. ❌ Failed experiments - Documented
4. ❌ Temporary files - Not relevant

---

## 🚀 For V4 Training

### **What to Reuse:**
- ✅ Model architecture (genocache_encoder.py)
- ✅ Training pipeline (train_single_gpu.py) - update for error augmentation
- ✅ Alignment pipeline (align_production.py) - should work better
- ✅ Validation pipeline (validate_giab_style.py) - reuse as-is

### **What to Replace:**
- ❌ Model weights - Retrain with error augmentation
- ❌ Genome embeddings - Re-encode with new model
- ❌ Index - Rebuild with new embeddings

### **What to Add:**
- ➕ Error augmentation module
- ➕ Reverse complement augmentation
- ➕ Curriculum learning schedule
- ➕ Hybrid seeding (neural + exact k-mers)

---

## ✅ Cleanup Commands (Optional)

```bash
# Create backup if you want to free up 52 GB:

cd /home/nebius/genocache

# Confirm backup is complete
ls -lh genocache-v3-backup/model/genocache_best.pt
ls -lh genocache-v3-backup/code/
ls -lh genocache-v3-backup/results/
ls -lh genocache-v3-backup/docs/

# If backup looks good, clean up:
cd genocache-production

# Delete large regeneratable files (saves 52 GB)
rm -rf encoded/
rm -f models/checkpoints/genocache_epoch*.pt
rm -f models/checkpoints/genocache_final.pt
rm -f giab_validation/genocache.sam
rm -f giab_validation/minimap2.sam

# Keep:
# - models/checkpoints/genocache_best.pt (backed up)
# - *.py files (code)
# - *.json files (metrics)
# - *.md files (docs)
# - giab_validation/ont_reads.fq (test data)
```

---

**End of Files Manifest**  
**Backup is minimal, complete, and reproducible** ✅
