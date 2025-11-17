#!/bin/bash
# Quick start script for GenoCache-Align

echo "🚀 GenoCache-Align Quick Start"
echo ""

# Check if conda env exists
if conda env list | grep -q neuralign; then
    echo "✅ Found neuralign environment"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate neuralign
else
    echo "Creating neuralign environment..."
    conda create -n neuralign python=3.10 -y
    conda activate neuralign
    pip install -r requirements.txt --break-system-packages
fi

echo ""
echo "Running validation test..."
python scripts/test_clustered_voting.py \
  --checkpoint checkpoints/improved_cnn_best.pt \
  --index indexes/index_ivfpq.faiss \
  --positions data/ref_positions_20251111_163637.npy \
  --fasta ../references/GCF_000001405.40_GRCh38.p14_chr22.fna \
  --chrom NC_000022.11 \
  --num-reads 100

echo ""
echo "✅ If you see ~99% accuracy, everything is working!"
