#!/bin/bash
# Project Organization Script
# Run this to clean up and prepare for download

echo "🧹 Organizing GenoCache-Align Project..."

# Create clean directory structure
mkdir -p organized/{scripts,indexes,data,checkpoints,outputs,docs}

# Copy working scripts only
echo "📁 Copying scripts..."
cp scripts/align_pipeline_final.py organized/scripts/
cp scripts/align_batched_fast.py organized/scripts/
cp scripts/align_batched_fp16.py organized/scripts/
cp scripts/test_clustered_voting.py organized/scripts/
cp scripts/test_clustered_with_errors.py organized/scripts/
cp scripts/debug_alignment.py organized/scripts/
cp scripts/improved_cnn.py organized/scripts/
cp scripts/genomic_utils.py organized/scripts/

# Copy essential data
echo "📊 Copying indexes and data..."
cp indexes/index_ivfpq.faiss organized/indexes/
cp data/reference_encodings/ref_positions_20251111_163637.npy organized/data/

# Copy checkpoint
echo "🧠 Copying model checkpoint..."
cp /home/nebius/work/genocache_checkpoints/improved_cnn_best.pt organized/checkpoints/

# Copy WFA tools
echo "🔧 Copying WFA tools..."
cp wfa_align organized/
cp wfa_align.cpp organized/

# Copy documentation
echo "📚 Copying documentation..."
cp README_FINAL.md organized/docs/README.md
cp RESULTS_SUMMARY.md organized/docs/
cp CONTINUE_NEXT_SESSION.md organized/docs/

# Copy example outputs
echo "📄 Copying example outputs..."
cp final_alignment_1k.sam organized/outputs/example_output.sam

# Create requirements file
cat > organized/requirements.txt << 'REQEOF'
torch>=2.0.0
numpy>=1.24.0
faiss-gpu>=1.7.4
biopython>=1.81
pysam>=0.21.0
REQEOF

# Create quick start script
cat > organized/quick_start.sh << 'QSEOF'
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
QSEOF

chmod +x organized/quick_start.sh

# Create tarball
echo "📦 Creating archive..."
cd organized
tar -czf ../genocache_v1.0.tar.gz *
cd ..

echo ""
echo "✅ Done! Files organized in: organized/"
echo "📦 Download: genocache_v1.0.tar.gz (~150MB)"
echo ""
echo "To download to local machine:"
echo "  scp nebius@<ip>:~/genocache/genocache_v2/genocache_v1.0.tar.gz ."
