#!/bin/bash
#
# Download GenoCache V4 artifacts from S3 for 8-GPU training
# Run this script in your NEW 8-GPU environment
#

set -e

# Configuration
S3_BUCKET="YOUR_BUCKET_NAME"  # ⚠️ CHANGE THIS!
S3_PREFIX="genocache-v4"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║     GenoCache V4 - Download from S3 for 8-GPU Training             ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "S3 Location: s3://${S3_BUCKET}/${S3_PREFIX}/"
echo ""

# Check AWS CLI
if ! command -v aws &> /dev/null; then
    echo "❌ AWS CLI not found! Install with: pip install awscli"
    exit 1
fi

# Test S3 access
echo "Testing S3 access..."
if ! aws s3 ls s3://${S3_BUCKET}/${S3_PREFIX}/ > /dev/null 2>&1; then
    echo "❌ Cannot access S3 bucket: ${S3_BUCKET}"
    echo "   Run: aws configure"
    exit 1
fi
echo "✅ S3 access confirmed"
echo ""

# Create directory structure
echo "Creating directory structure..."
mkdir -p {models,scripts,validation,logs,data}
echo "✅ Directories created"
echo ""

# Download setup guide first
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Downloading setup guide..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
aws s3 cp s3://${S3_BUCKET}/${S3_PREFIX}/SETUP_GUIDE.md ./SETUP_GUIDE.md
echo "✅ Setup guide downloaded"
echo ""

# Download genome
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Downloading genome reference (3.1 GB)..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
aws s3 cp s3://${S3_BUCKET}/${S3_PREFIX}/genome/GRCh38.fa ./GRCh38.fa
aws s3 cp s3://${S3_BUCKET}/${S3_PREFIX}/genome/GRCh38.fa.fai ./GRCh38.fa.fai
echo "✅ Genome downloaded"
echo ""

# Download scripts
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Downloading training scripts..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
aws s3 sync s3://${S3_BUCKET}/${S3_PREFIX}/scripts/ ./scripts/
chmod +x ./scripts/*.sh 2>/dev/null || true
echo "✅ Scripts downloaded"
echo ""

# Download model (optional, for reference or warm start)
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Downloading current best model (for reference)..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
aws s3 sync s3://${S3_BUCKET}/${S3_PREFIX}/models/ ./models/ --exclude "*" --include "*.pt"
echo "✅ Model(s) downloaded"
echo ""

# Download validation data (optional)
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Downloading validation data..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
aws s3 sync s3://${S3_BUCKET}/${S3_PREFIX}/validation/ ./validation/
echo "✅ Validation data downloaded"
echo ""

# Summary
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    DOWNLOAD COMPLETE! ✅                            ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Files downloaded:"
echo "  ✅ GRCh38.fa (genome reference)"
echo "  ✅ scripts/ (training code)"
echo "  ✅ models/ (reference model)"
echo "  ✅ validation/ (test data)"
echo "  ✅ SETUP_GUIDE.md (instructions)"
echo ""
echo "Next steps:"
echo "  1. Read the setup guide:"
echo "     cat SETUP_GUIDE.md"
echo ""
echo "  2. Setup environment:"
echo "     conda create -n genocache python=3.11 -y"
echo "     conda activate genocache"
echo "     pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121"
echo "     pip install biopython numpy tqdm faiss-gpu pandas matplotlib seaborn"
echo ""
echo "  3. Test PyTorch:"
echo "     python3 -c 'import torch; print(f\"GPUs: {torch.cuda.device_count()}\")'"
echo ""
echo "  4. Start training:"
echo "     torchrun --nproc_per_node=8 scripts/train_8gpu.py \\"
echo "       --genome GRCh38.fa \\"
echo "       --batch-size 1024 \\"
echo "       --epochs 50 \\"
echo "       --output-dir ./models/8gpu"
echo ""
echo "🚀 Ready to train with batch=8192 (8 GPUs × 1024)!"
echo ""
