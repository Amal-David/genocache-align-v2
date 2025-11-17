# GenoCache Version 6.1 - Full Genome NAL Training

**Created:** November 16, 2025  
**Status:** Ready for full genome training

## What's New in v6.1

### Based on Chr22 Success (v4.1)
- Chr22 focused training: 76.5% accuracy (beats minimap2!)
- Proved that more training + curriculum learning works

### Scaling to Full Genome
- **Training:** ALL chromosomes (GRCh38, ~3 Gbp)
- **Batches:** 16,000 (2x more than Chr22's 8000, conservative for GPU)
- **Curriculum:** Same proven 5% → 15% error progression
- **Architecture:** Same CNN-based NAL encoder
- **Seed length:** 512bp (not using 256/512 variants)
- **Multi-GPU:** Auto-detects and uses all available GPUs
- **Expected time:** 3-5 hours of training

## Directory Structure

```
GenoCache_Version6.1/
├── training/           # Full genome training scripts
│   ├── train_full_genome.py       # Main training (24K batches)
│   └── encoder_nal.py             # Model architecture
├── indexing/           # Full genome index building
│   └── build_full_genome_index.py # Build complete index
├── alignment/          # Alignment pipeline
│   ├── seeding_nal.py             # Seed extraction
│   ├── chaining_nal.py            # Anchor chaining
│   └── test_alignment.py          # Testing
├── models/             # Trained models (created during training)
├── indexes/            # FAISS indexes (created during indexing)
├── results/            # Test results
├── data/               # Test datasets
└── documentation/      # All documentation

## Key Parameters

### Training
- Chromosomes: All 24 (Chr1-22, X, Y)
- Batches: 16,000 (conservative for GPU memory)
- Batch size: 32 (auto-adjusts if GPU memory < 40GB)
- Curriculum: 5% → 8% → 12% → 15% error rate
- Seed length: 512bp (single length, proven in Chr22)
- Embedding dim: 128
- Error types: Substitution, Insertion, Deletion (independent)
- Multi-GPU: Auto-detects and parallelizes across all GPUs

### Expected Results
- Target accuracy: 75-85% (based on Chr22 76.5%)
- Training time: ~4-6 hours
- Index building: ~2-3 hours
- Final index size: ~20-30 GB

## How to Use

### 1. Train Model
```bash
cd training/
python3 train_full_genome.py
# Takes ~4-6 hours
# Creates: models/genocache_nal_fullgenome.pt
```

### 2. Build Index
```bash
cd ../indexing/
python3 build_full_genome_index.py
# Takes ~2-3 hours
# Creates: indexes/genocache_nal_fullgenome.{index,positions.npz}
```

### 3. Test Alignment
```bash
cd ../alignment/
python3 test_alignment.py
# Tests on synthetic/real reads
```

## Differences from v4.1 (Chr22)

| Feature | v4.1 (Chr22) | v6.1 (Full Genome) |
|---------|--------------|-------------------|
| Chromosomes | 1 (Chr22) | 24 (All) |
| Genome size | 50.8 Mbp | ~3 Gbp (60x larger) |
| Batches | 8,000 | 16,000 (2x more) |
| Training time | ~1.5 hours | ~3-5 hours (2x) |
| Index size | 1.1 GB | ~20-30 GB (20x) |
| Expected accuracy | 76.5% | 75-85% |

## Requirements

- Python 3.12+
- PyTorch with CUDA
- FAISS-GPU
- ~40 GB free disk space
- NVIDIA GPU (H100 recommended)
- GRCh38 reference genome

## Notes

- Uses same proven curriculum learning from Chr22
- Same error injection (Sub/Ins/Del) as NeuralAligner paper
- Scales the successful Chr22 approach to full genome
- Training samples from all chromosomes uniformly
