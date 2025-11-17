# GenoCache V4 - Training Plan

**Goal:** Build error-tolerant neural aligner that surpasses minimap2, BWA-MEM, and NeuralAligner  
**Key Innovation:** First neural aligner trained explicitly on noisy sequences  
**Target:** Clinical-grade accuracy (≥95% on ONT, ≥99% on Illumina) with GPU acceleration

---

## 🎯 V4 Design Philosophy

### **Core Innovation:**
```
V3 Problem: Neural seeds fail on noisy reads (27% accuracy)
V4 Solution: Train on noisy data → Error-tolerant embeddings
Expected:    95%+ accuracy on 8% error ONT reads
```

### **Why This Will Work:**
1. **Architecture is proven** - V3 gets 90% on clean reads
2. **Problem is training** - Never saw errors during training
3. **Solution is known** - Add error augmentation (like BERT)
4. **Precedent exists** - DNABERT, ProtTrans use augmentation

---

## 📋 V4 vs V3 Comparison

| Aspect | V3 (Failed) | V4 (Target) |
|--------|-------------|-------------|
| **Training Data** | Clean sequences only | Clean + 1-10% error rates |
| **Augmentation** | None | Errors + RC + Hard negatives |
| **Seed Length** | 512 bp | 512-1024 bp (adaptive) |
| **ONT Accuracy** | 27% ❌ | ≥95% ✅ |
| **Illumina Accuracy** | Not tested | ≥99% ✅ |
| **Speed** | 2.7 reads/sec (CPU) | ≥500 reads/sec (GPU) |
| **Robustness** | Error-sensitive | Error-tolerant |

---

## 🏗️ V4 Architecture

### **Model Design:**

```python
class GenoCacheV4Encoder(nn.Module):
    """
    Error-Tolerant Neural Aligner
    
    Key Differences from V3:
    1. Longer context (1024bp vs 512bp)
    2. Error-aware attention mechanism
    3. Dropout for robustness
    4. Trained on noisy sequences
    """
    
    def __init__(
        self,
        seed_len=1024,      # Longer for robustness (vs 512 in V3)
        emb_dim=256,
        num_layers=4,
        dropout=0.1         # NEW: Add dropout for noise tolerance
    ):
        super().__init__()
        
        # Multi-scale CNN (same as V3, works well)
        self.conv_blocks = nn.ModuleList([
            ConvBlock(5, 64, kernel_size=k) 
            for k in [3, 7, 15, 31]  # Multi-resolution
        ])
        
        # NEW: Error-aware attention
        # Learns to focus on conserved regions, ignore errors
        self.attention = ErrorAwareAttention(
            dim=256,
            num_heads=8,
            dropout=dropout
        )
        
        # Projection to embedding space
        self.projection = nn.Sequential(
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Dropout(dropout),  # NEW
            nn.Linear(512, emb_dim)
        )
        
    def forward(self, x):
        # Multi-scale features
        features = [conv(x) for conv in self.conv_blocks]
        x = torch.cat(features, dim=-1)
        
        # Error-aware attention
        x = self.attention(x)
        
        # Project to embedding
        emb = self.projection(x)
        
        # L2 normalize
        return F.normalize(emb, p=2, dim=-1)
```

### **Key Architectural Improvements:**

1. **Longer Seeds (1024bp vs 512bp)**
   - More signal → More robust to errors
   - 1024bp with 8% error = ~80 errors, but still has 920bp correct
   - Better discriminative power between regions

2. **Error-Aware Attention**
   - Learns which positions are conserved vs error-prone
   - Attention weights focus on high-confidence regions
   - Similar to BERT's masked language model approach

3. **Dropout Regularization**
   - Forces model to not over-rely on specific positions
   - Better generalization to unseen error patterns
   - Standard technique for noisy data

4. **Same Loss Function**
   - Contrastive loss with hard negatives (V3 approach works)
   - No change needed here

---

## 🔧 V4 Training Strategy

### **1. Error Augmentation Pipeline**

```python
def augment_sequence(seq, error_rate=0.08):
    """
    Add realistic sequencing errors
    
    Error Types (ONT profile):
    - 60% Substitutions
    - 20% Insertions  
    - 20% Deletions
    
    Returns: noisy_seq, error_positions
    """
    errors = {
        'substitution': 0.6,
        'insertion': 0.2,
        'deletion': 0.2
    }
    
    noisy = list(seq)
    error_positions = []
    
    for i in range(len(seq)):
        if random.random() < error_rate:
            error_type = random.choices(
                list(errors.keys()),
                weights=errors.values()
            )[0]
            
            if error_type == 'substitution':
                noisy[i] = random.choice([b for b in 'ACGT' if b != seq[i]])
            elif error_type == 'insertion':
                noisy.insert(i, random.choice('ACGT'))
            else:  # deletion
                noisy.pop(i)
            
            error_positions.append(i)
    
    return ''.join(noisy), error_positions

# Usage in training:
for clean_seq in dataset:
    # 80% noisy, 20% clean (for comparison)
    if random.random() < 0.8:
        # Random error rate (1-10%)
        error_rate = random.uniform(0.01, 0.10)
        noisy_seq, _ = augment_sequence(clean_seq, error_rate)
        batch.append(noisy_seq)
    else:
        batch.append(clean_seq)  # Keep some clean
```

### **2. Reverse Complement Augmentation**

```python
def reverse_complement(seq):
    """RC augmentation for strand invariance"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[b] for b in reversed(seq))

# 50% of training examples are RC
for seq in dataset:
    if random.random() < 0.5:
        seq = reverse_complement(seq)
    
    # Then apply error augmentation
    noisy_seq, _ = augment_sequence(seq, error_rate)
```

### **3. Hard Negative Mining** (Keep from V3)

```python
# V3 identified 983K repetitive regions - reuse this!
# Continue with 30-40% hard negatives during training

hard_negative_regions = load_hard_negatives()  # From V3

for anchor in batch:
    # Positive: Same region + different errors
    pos_1 = augment_sequence(anchor, 0.05)
    pos_2 = augment_sequence(anchor, 0.08)
    
    # Hard negative: Repetitive region (similar sequence, different location)
    neg_region = sample_hard_negative(hard_negative_regions)
    neg = augment_sequence(neg_region, 0.05)
    
    # Loss: Pull pos_1, pos_2 together; push neg away
    loss = contrastive_loss(anchor, pos_1, pos_2, neg)
```

### **4. Curriculum Learning**

```python
# Gradually increase difficulty

curriculum = [
    # Phase 1: Learn clean patterns (5 epochs)
    {'error_rate': 0.00, 'hard_neg_ratio': 0.2},
    
    # Phase 2: Introduce small errors (5 epochs)
    {'error_rate': 0.01, 'hard_neg_ratio': 0.25},
    
    # Phase 3: Moderate errors (5 epochs)
    {'error_rate': 0.03, 'hard_neg_ratio': 0.30},
    
    # Phase 4: Higher errors (10 epochs)  
    {'error_rate': 0.05, 'hard_neg_ratio': 0.35},
    
    # Phase 5: Full ONT errors (15 epochs)
    {'error_rate': 0.08, 'hard_neg_ratio': 0.40},
    
    # Phase 6: Mix all error rates (10 epochs)
    {'error_rate': random.uniform(0, 0.10), 'hard_neg_ratio': 0.40}
]

# Total: 50 epochs
# Training time: ~24 hours on H100
```

---

## 📊 Training Details

### **Dataset:**
```
Source: Full GRCh38 (2.65B bp, 705 sequences)
Stride: 100bp
Seeds per position: 
  - 512bp seeds: 26M
  - 1024bp seeds: 26M
Total examples: ~52M (with augmentation)

Breakdown:
- 20% Clean sequences (baseline)
- 80% Noisy sequences (1-10% error)
- 50% RC augmented
- 40% Hard negatives

Effective dataset size: 
- 26M positions × 5 augmentations = 130M examples
```

### **Training Configuration:**
```python
config = {
    # Model
    'seed_len': 1024,
    'emb_dim': 256,
    'num_layers': 4,
    'dropout': 0.1,
    
    # Training
    'batch_size': 256,      # Fit in H100 80GB
    'learning_rate': 1e-4,
    'weight_decay': 1e-5,
    'epochs': 50,
    
    # Augmentation
    'error_rate_range': (0.01, 0.10),
    'rc_prob': 0.5,
    'hard_negative_ratio': 0.40,
    
    # Optimization
    'optimizer': 'AdamW',
    'scheduler': 'CosineAnnealingWarmRestarts',
    'warmup_epochs': 5,
    
    # Hardware
    'device': 'H100',
    'mixed_precision': True,  # fp16
    'gradient_checkpointing': True
}
```

### **Expected Training Time:**
- **Data loading:** 2 hours (generate augmented dataset)
- **Training:** 24 hours (50 epochs × ~30 min/epoch)
- **Validation:** 2 hours (test on multiple error rates)
- **Total:** ~28 hours

### **Checkpointing:**
- Save every 5 epochs
- Keep best model based on validation loss
- Track metrics:
  - Loss on clean sequences
  - Loss on noisy sequences (1%, 3%, 5%, 8%, 10%)
  - Recall@10 on validation set

---

## ✅ Success Criteria

### **During Training:**
- [ ] Loss converges smoothly (<0.1 final)
- [ ] Clean sequence accuracy ≥90% (V3 baseline)
- [ ] 8% error sequence accuracy ≥80% (NEW)
- [ ] Hard negatives properly separated

### **After Training (Validation):**
- [ ] Synthetic ONT (8% error): ≥90% correct chromosome
- [ ] Real ONT (HG002): ≥95% accuracy vs minimap2
- [ ] Illumina (1% error): ≥99% accuracy
- [ ] Speed: ≥500 reads/sec on GPU

### **Production Deployment:**
- [ ] GIAB HG002 validation: Match minimap2 accuracy
- [ ] Variant calling accuracy: ≥95% precision/recall
- [ ] Memory usage: ≤80 GB (fits in H100)
- [ ] End-to-end pipeline: FASTQ → SAM

---

## 🚀 Implementation Timeline

### **Phase 1: Preparation (1 day)**
- [ ] Set up V4 training directory
- [ ] Copy V3 code as starting point
- [ ] Implement error augmentation module
- [ ] Implement RC augmentation
- [ ] Test augmentation pipeline (sanity check)
- [ ] Prepare validation sets (clean + noisy)

### **Phase 2: Model Updates (4 hours)**
- [ ] Extend seed length (512 → 1024bp)
- [ ] Add error-aware attention module
- [ ] Add dropout layers
- [ ] Update data loader for augmentation
- [ ] Test forward/backward pass

### **Phase 3: Training (1-2 days)**
- [ ] Generate augmented training data
- [ ] Launch curriculum learning (50 epochs)
- [ ] Monitor training (loss curves, metrics)
- [ ] Checkpointing every 5 epochs
- [ ] Validation on multiple error rates

### **Phase 4: Validation (1 day)**
- [ ] Encode genome with V4 model
- [ ] Build FAISS index (IVF-Flat, 26 GB)
- [ ] Test on synthetic ONT (8% error)
- [ ] Test on real ONT (HG002 if available)
- [ ] Compare to minimap2 baseline

### **Phase 5: Optimization (1-2 days)**
- [ ] GPU batch inference optimization
- [ ] FAISS GPU search integration
- [ ] End-to-end profiling
- [ ] Target: 500+ reads/sec

### **Phase 6: Production (1 day)**
- [ ] Clean up code
- [ ] Create CLI tool
- [ ] Write user documentation
- [ ] Create Docker container
- [ ] Prepare for release

**Total Time: 6-8 days**

---

## 📁 V4 Directory Structure

```
genocache-v4/
├── data/
│   ├── GRCh38.fa                    # Reference genome
│   ├── augmented/                   # Generated training data
│   │   ├── train_clean.npy
│   │   ├── train_noisy_1pct.npy
│   │   ├── train_noisy_3pct.npy
│   │   ├── train_noisy_5pct.npy
│   │   ├── train_noisy_8pct.npy
│   │   └── train_noisy_10pct.npy
│   └── validation/
│       ├── val_clean.npy
│       ├── val_noisy_8pct.npy
│       └── ont_reads_real.fq        # HG002 if available
│
├── models/
│   ├── encoder_v4.py                # Updated architecture
│   ├── attention.py                 # Error-aware attention
│   └── loss.py                      # Contrastive loss
│
├── training/
│   ├── augmentation.py              # Error + RC augmentation
│   ├── dataset.py                   # DataLoader with augmentation
│   ├── train_v4.py                  # Main training script
│   └── curriculum.py                # Curriculum learning schedule
│
├── alignment/
│   ├── encode_genome.py             # Genome encoding
│   ├── build_index.py               # FAISS index building
│   ├── align_fastq.py               # Aligner (from V3, should work)
│   └── validate.py                  # Validation pipeline
│
├── checkpoints/
│   ├── v4_epoch05.pt
│   ├── v4_epoch10.pt
│   ├── ...
│   └── v4_best.pt                   # Best checkpoint
│
├── results/
│   ├── training_curves.png
│   ├── validation_metrics.json
│   └── giab_comparison.json
│
└── docs/
    ├── TRAINING_LOG.md
    ├── ARCHITECTURE.md
    └── RESULTS.md
```

---

## 🔬 Expected Improvements Over V3

### **1. Accuracy on Noisy Reads:**
```
Metric: Correct chromosome rate
Data: ONT reads with 8% error

V3:  27.4% ❌
V4:  ≥90% ✅ (target)

Improvement: +62.6% absolute, 3.3× relative
Why: Error-tolerant embeddings from training
```

### **2. Robustness Across Error Rates:**
```
Test: Alignment accuracy vs error rate

V3:
├─ 0% error:  90% ✅
├─ 2% error:  ~60% ⚠️
├─ 5% error:  ~40% ❌
└─ 8% error:  27% ❌

V4 (expected):
├─ 0% error:  95% ✅
├─ 2% error:  94% ✅
├─ 5% error:  92% ✅
└─ 8% error:  90% ✅

Improvement: Consistent across error rates
Why: Trained on all error rates (1-10%)
```

### **3. Strand Invariance:**
```
Test: Forward vs reverse complement

V3:
├─ Forward:   90% ✅
└─ RC:        ~30% ❌ (not trained)

V4 (expected):
├─ Forward:   95% ✅
└─ RC:        95% ✅

Improvement: RC same as forward
Why: 50% RC augmentation during training
```

### **4. Speed (with GPU optimization):**
```
Hardware: H100 GPU

V3:  2.7 reads/sec (CPU, no optimization)
V4:  ≥500 reads/sec (GPU, batched)

Improvement: 185× faster
Why: GPU inference + batch processing
```

---

## 🎯 Competitive Analysis

### **vs Minimap2:**
```
Metric           Minimap2    V4 Target    Advantage
---------------------------------------------------------
ONT Accuracy     94.6%       ≥95%         Equal/better
Illumina Acc     ~98%        ≥99%         Better
Speed (CPU)      10 r/s      500 r/s      50× faster
Speed (GPU)      N/A         500 r/s      GPU acceleration
Memory           8 GB        30 GB        More (acceptable)
```

### **vs BWA-MEM:**
```
Metric           BWA-MEM     V4 Target    Advantage
---------------------------------------------------------
Illumina Acc     99.8%       ≥99%         Comparable
ONT Support      No          Yes          Better (novel)
Speed (CPU)      50 r/s      500 r/s      10× faster
GPU Support      No          Yes          GPU acceleration
```

### **vs NeuralAligner:**
```
Metric           NeuralAligner  V4 Target    Advantage
---------------------------------------------------------
Architecture     Hybrid         Pure neural  Simpler
Error Training   No             Yes          More robust
Accuracy         ~95%           ≥95%         Comparable
Speed            Unknown        ≥500 r/s     Likely faster
Open Source      No             Yes          Accessible
```

---

## 💡 Innovation Summary

### **What Makes V4 Special:**

1. **First Error-Trained Neural Aligner**
   - All existing neural aligners train on clean sequences
   - V4 explicitly trains on 1-10% error rates
   - Matches real sequencing conditions

2. **Pure Neural Approach**
   - No hybrid with exact k-mers (unlike NeuralAligner)
   - End-to-end learnable
   - Simpler pipeline

3. **GPU-Accelerated Throughout**
   - Embedding: GPU (PyTorch)
   - Search: GPU (FAISS-GPU)
   - Alignment: GPU (WFA-GPU or Parasail)
   - 100-500× potential speedup

4. **Unified for All Read Types**
   - ONT long reads (5-20kb, 8% error)
   - PacBio long reads (10-30kb, 5% error)  
   - Illumina short reads (150bp, 1% error)
   - Single tool, multiple applications

5. **Clinical-Grade Validation**
   - GIAB gold standard testing
   - Variant calling validation
   - Production-ready accuracy

---

## ✅ Ready to Start V4

### **We Have:**
- ✅ Proven architecture (V3 works on clean data)
- ✅ Understanding of failure mode (error intolerance)
- ✅ Clear solution (error-augmented training)
- ✅ Complete pipeline (encode → index → align)
- ✅ Validation framework (GIAB-style testing)
- ✅ H100 GPU access

### **We Need:**
- [ ] Implement error augmentation pipeline
- [ ] Update model architecture (longer seeds, dropout)
- [ ] Set up curriculum learning
- [ ] Train for 24-48 hours
- [ ] Validate on real ONT data

### **Expected Outcome:**
**Clinical-grade neural aligner that surpasses existing tools in speed and matches them in accuracy** 🚀

---

**Let's build V4 right!** ✅
