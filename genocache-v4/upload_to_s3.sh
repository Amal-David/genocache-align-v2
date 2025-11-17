#!/bin/bash
#
# Upload GenoCache V4 artifacts to S3 for 8-GPU training
#

set -e

# Configuration
S3_BUCKET="YOUR_BUCKET_NAME"  # ⚠️ CHANGE THIS!
S3_PREFIX="genocache-v4"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║        Uploading GenoCache V4 Artifacts to S3                      ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "S3 Bucket: s3://${S3_BUCKET}/${S3_PREFIX}/"
echo ""

# Check AWS CLI
if ! command -v aws &> /dev/null; then
    echo "❌ AWS CLI not found! Install with: pip install awscli"
    exit 1
fi

# Test S3 access
echo "Testing S3 access..."
if ! aws s3 ls s3://${S3_BUCKET}/ > /dev/null 2>&1; then
    echo "❌ Cannot access S3 bucket: ${S3_BUCKET}"
    echo "   Run: aws configure"
    exit 1
fi
echo "✅ S3 access confirmed"
echo ""

# ============================================================================
# PART 1: Genome Reference (Essential)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PART 1: Uploading Genome Reference"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "/home/nebius/genocache/GRCh38.fa" ]; then
    echo "Uploading GRCh38.fa (3.1 GB)..."
    aws s3 cp /home/nebius/genocache/GRCh38.fa \
        s3://${S3_BUCKET}/${S3_PREFIX}/genome/GRCh38.fa \
        --storage-class INTELLIGENT_TIERING
    echo "✅ GRCh38.fa uploaded"
else
    echo "⚠️  GRCh38.fa not found in /home/nebius/genocache/"
fi

if [ -f "/home/nebius/genocache/GRCh38.fa.fai" ]; then
    echo "Uploading GRCh38.fa.fai..."
    aws s3 cp /home/nebius/genocache/GRCh38.fa.fai \
        s3://${S3_BUCKET}/${S3_PREFIX}/genome/GRCh38.fa.fai
    echo "✅ GRCh38.fa.fai uploaded"
fi

echo ""

# ============================================================================
# PART 2: Best Model Checkpoint (Essential)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PART 2: Uploading Model Checkpoints"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -d "/home/nebius/genocache/genocache-v4/models/checkpoints" ]; then
    echo "Uploading all model checkpoints..."
    aws s3 sync /home/nebius/genocache/genocache-v4/models/checkpoints/ \
        s3://${S3_BUCKET}/${S3_PREFIX}/models/ \
        --exclude "*" \
        --include "*.pt"
    echo "✅ Model checkpoints uploaded"
else
    echo "⚠️  Model directory not found"
fi

echo ""

# ============================================================================
# PART 3: Training Scripts (Essential)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PART 3: Uploading Training Scripts"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Create temporary scripts directory
SCRIPTS_DIR="/tmp/genocache_scripts_$$"
mkdir -p ${SCRIPTS_DIR}

# Copy all Python scripts
if [ -d "/home/nebius/genocache/genocache-v4" ]; then
    cp /home/nebius/genocache/genocache-v4/model.py ${SCRIPTS_DIR}/ 2>/dev/null || true
    cp /home/nebius/genocache/genocache-v4/dataset.py ${SCRIPTS_DIR}/ 2>/dev/null || true
    cp /home/nebius/genocache/genocache-v4/train.py ${SCRIPTS_DIR}/ 2>/dev/null || true
    cp /home/nebius/genocache/genocache-v4/train_8gpu.py ${SCRIPTS_DIR}/ 2>/dev/null || true
    cp /home/nebius/genocache/genocache-v4/validate.py ${SCRIPTS_DIR}/ 2>/dev/null || true
    cp /home/nebius/genocache/genocache-v4/validation/*.py ${SCRIPTS_DIR}/ 2>/dev/null || true
fi

echo "Uploading training scripts..."
aws s3 sync ${SCRIPTS_DIR}/ s3://${S3_BUCKET}/${S3_PREFIX}/scripts/ \
    --exclude "*" \
    --include "*.py"

# Upload setup guide
if [ -f "/home/nebius/genocache/genocache-v4/SETUP_NEW_TRAINING_8GPU.md" ]; then
    aws s3 cp /home/nebius/genocache/genocache-v4/SETUP_NEW_TRAINING_8GPU.md \
        s3://${S3_BUCKET}/${S3_PREFIX}/SETUP_GUIDE.md
    echo "✅ Setup guide uploaded"
fi

rm -rf ${SCRIPTS_DIR}
echo "✅ Training scripts uploaded"
echo ""

# ============================================================================
# PART 4: Validation Data (Optional but Useful)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PART 4: Uploading Validation Data"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -d "/home/nebius/genocache/genocache-v4/validation" ]; then
    echo "Uploading validation data..."
    aws s3 sync /home/nebius/genocache/genocache-v4/validation/ \
        s3://${S3_BUCKET}/${S3_PREFIX}/validation/ \
        --exclude "*.index" \
        --exclude "*.pkl"
    echo "✅ Validation data uploaded"
else
    echo "⚠️  Validation directory not found"
fi

echo ""

# ============================================================================
# PART 5: Documentation and Logs (Optional)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PART 5: Uploading Documentation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Upload key documentation files
for doc in FINAL_RESULTS.md NEXT_PHASE_PLAN.md README_BENCHMARK.md; do
    if [ -f "/home/nebius/genocache/genocache-v4/${doc}" ]; then
        aws s3 cp "/home/nebius/genocache/genocache-v4/${doc}" \
            "s3://${S3_BUCKET}/${S3_PREFIX}/docs/${doc}"
        echo "✅ ${doc} uploaded"
    fi
done

echo ""

# ============================================================================
# SUMMARY
# ============================================================================
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    UPLOAD COMPLETE! ✅                              ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Files uploaded to: s3://${S3_BUCKET}/${S3_PREFIX}/"
echo ""
echo "Structure:"
echo "  s3://${S3_BUCKET}/${S3_PREFIX}/"
echo "  ├── genome/"
echo "  │   ├── GRCh38.fa"
echo "  │   └── GRCh38.fa.fai"
echo "  ├── models/"
echo "  │   └── *.pt"
echo "  ├── scripts/"
echo "  │   ├── model.py"
echo "  │   ├── dataset.py"
echo "  │   ├── train_8gpu.py"
echo "  │   └── validate.py"
echo "  ├── validation/"
echo "  │   └── *.json, *.fa"
echo "  ├── docs/"
echo "  │   └── *.md"
echo "  └── SETUP_GUIDE.md"
echo ""
echo "Next steps:"
echo "  1. In new 8-GPU environment, run:"
echo "     aws s3 cp s3://${S3_BUCKET}/${S3_PREFIX}/SETUP_GUIDE.md ."
echo "     cat SETUP_GUIDE.md"
echo ""
echo "  2. Follow the guide to download files and start training!"
echo ""
echo "You can share this with new Droid by saying:"
echo "  'Download from s3://${S3_BUCKET}/${S3_PREFIX}/ and follow SETUP_GUIDE.md'"
echo ""
