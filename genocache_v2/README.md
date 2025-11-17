# GenoCache-Align V2 - Clean Architecture

## 🎯 Vision
Two-tier genomic aligner:
- **Tier 1:** Fast CNN for ~70-80% of reads (unique mappers)
- **Tier 2:** Hyena rescue for repetitive/homologous regions

## 📁 Directory Structure

```
genocache_v2/
├── models/
│   ├── cnn_fast/              # Fast CNN tier (ImprovedCNN)
│   │   ├── architecture.py
│   │   ├── train.py
│   │   └── model_manifest.json
│   ├── hyena_rescue/          # Hyena rescue tier (future)
│   └── checkpoints/           # Saved model weights
├── data/
│   ├── reference/             # Reference genome embeddings
│   ├── queries/               # Query read embeddings
│   └── manifests/             # Provenance tracking
├── indexes/
│   ├── flat/                  # Uncompressed FAISS (validation)
│   ├── ivf/                   # IVF partitioned
│   └── compressed/            # IVFPQ/OPQ (production)
├── scripts/
│   ├── encode_reference.py   # Reference → embeddings
│   ├── encode_queries.py     # Reads → embeddings
│   └── build_index.py        # Embeddings → FAISS
└── logs/                      # Training/encoding logs
```

## 🔑 Key Principles
1. **One model, one checkpoint** - No architecture confusion
2. **Manifest everything** - Track which checkpoint created which embeddings
3. **Version control** - Date-stamped outputs
4. **GPU-optimized** - Built for 8x GPU scaling
