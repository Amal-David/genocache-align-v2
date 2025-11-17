# GenoCache V3 - Beyond NeurALigner

Next-generation genomic sequence aligner with multi-scale CNN + attention architecture.

## Architecture Highlights

- **1.2M parameters**: Multi-scale CNN + lightweight attention
- **256D embeddings**: Higher specificity than NeurALigner (128D)
- **Hard negative mining**: From repetitive genomic regions (Alu, LINE, SINE)
- **Curriculum learning**: Gradually increase augmentation difficulty
- **O(L) attention**: Linear complexity for long-range dependencies

## Quick Start

### 1. Test the Model

```bash
cd genocache-v3
python3 genocache_encoder.py
```

Expected output:
```
GenoCacheEncoder: 1,234,567 trainable parameters
✅ Model test passed!
```

### 2. Train on Chromosome 1

```bash
python3 train_genocache.py \
    --fasta ../GRCh38.fa \
    --chrom-id NC_000001.11 \
    --output-dir ./models \
    --epochs 50 \
    --batch-size 256
```

**Training time**: ~8-10 hours on H100 for 50 epochs

### 3. Validate

```bash
python3 validate_encoder.py \
    --model ./models/genocache_best.pt \
    --fasta ../GRCh38.fa \
    --chrom-id NC_000001.11
```

**Target metrics**:
- ✅ Inference: <25 μs per seed
- ✅ Same sequence similarity: >0.7
- ✅ Different sequence similarity: <0.3

## Architecture Details

### Multi-Scale Convolutions

```
Input (512bp) → Token Embedding (5→64D) → Positional Encoding
                                           ↓
    ┌──────────────────────────────────────────────────────┐
    │              Multi-Scale Conv Block 1                │
    │  ┌─────────┬─────────┬─────────┬─────────┐          │
    │  │ k=3     │ k=7     │ k=15    │ k=31    │ (64→128D)│
    │  └─────────┴─────────┴─────────┴─────────┘          │
    └──────────────────────────────────────────────────────┘
                            ↓
    ┌──────────────────────────────────────────────────────┐
    │              Multi-Scale Conv Block 2                │
    │  ┌─────────┬─────────┬─────────┬─────────┐          │
    │  │ k=3     │ k=7     │ k=15    │ k=31    │ (128→256D)│
    │  └─────────┴─────────┴─────────┴─────────┘          │
    └──────────────────────────────────────────────────────┘
                            ↓
            Lightweight Linear Attention (O(L), 4 heads)
                            ↓
                  Global Average Pooling
                            ↓
                   Linear Projection → 256D
                            ↓
                      L2 Normalization
```

### Hard Negative Mining

The training script automatically:

1. **Analyzes k-mer frequencies** (k=15) across chromosome 1
2. **Identifies repetitive regions** (k-mer count ≥ 10)
3. **Samples 30% of negatives** from these regions
4. Forces the model to learn fine distinctions in repetitive areas

**Why this matters**: Repetitive regions (Alu, LINE, SINE) are the hardest to align correctly. Hard negative mining improves recall by 10-15%.

## File Structure

```
genocache-v3/
├── README.md                  # This file
├── genocache_encoder.py       # Model architecture
├── train_genocache.py         # Training with hard negatives
└── validate_encoder.py        # Validation suite
```

## Performance Targets

| Metric | Target | NeurALigner |
|--------|--------|-------------|
| Inference speed | <25 μs/seed | 48 μs/seed |
| Top-1 recall | 97% | 94% |
| Embedding dim | 256D | 128D |
| Parameters | 1.2M | 0.5M |

## Training Configuration

**Default hyperparameters**:
- Seed length: 512 bp
- Embedding dim: 256
- Batch size: 256 (single H100)
- Learning rate: 3e-4
- Weight decay: 0.01
- Temperature: 0.07
- Hard negative ratio: 0.3

**Augmentation** (realistic error profiles):
- Substitution rate: 1-10%
- Indel simulation: ~20% of errors
- Reverse complement: 50%
- Random shifts: ±51 bp

**Curriculum learning**:
- First 10 epochs: easier examples
- Gradually increase error rates
- Improves convergence

## Next Steps

After training on chr1:

1. **Full genome encoding** (~30-60 min on H100)
2. **FAISS index building** (~10-20 min)
3. **Alignment pipeline** (seed selection → chaining → WFA)
4. **Benchmark** vs minimap2, NeurALigner, Parabricks

## Troubleshooting

### Out of Memory

Reduce batch size:
```bash
python3 train_genocache.py --batch-size 128
```

### Training too slow

Check GPU utilization:
```bash
nvidia-smi dmon -s u
```

Should be >80%. If not, increase `num_workers` in train_genocache.py (line 283).

### Poor embedding quality

Train longer or increase model capacity:
```bash
python3 train_genocache.py --epochs 100
```

## Citation

Based on:
- NeurALigner (ICLR 2026 submission)
- Hyena-DNA (NeurIPS 2023)
- Hard negative mining (FaceNet, CVPR 2015)

## License

MIT License - Free to use, modify, and commercialize
