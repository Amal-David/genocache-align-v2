#!/bin/bash
# Comprehensive backup of EVERYTHING to S3

# S3 Configuration
export AWS_ACCESS_KEY_ID="NAKIADXF5GT3EIV7XGCC"
export AWS_SECRET_ACCESS_KEY="b2hCHc4rpq+/HWqcF2tB61ATRRRlA9Uf34qEUICt"
export AWS_DEFAULT_REGION="us-east-1"
BUCKET="hackathon-team-fabric3-5"
PREFIX="genocache_full_backup"

echo "================================================================================"
echo "COMPREHENSIVE BACKUP TO S3 - ALL FILES"
echo "================================================================================"
echo "Bucket: s3://$BUCKET/$PREFIX/"
echo ""

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
echo "Backup timestamp: $TIMESTAMP"
echo ""

# 1. Backup ALL model checkpoints (415 files)
echo "📦 [1/10] Backing up ALL model checkpoints (.pt files)..."
aws s3 sync . s3://$BUCKET/$PREFIX/models/ \
    --exclude "*" \
    --include "*.pt" \
    --no-progress &
PID1=$!

# 2. Backup H100 vectors (148 GB)
echo "📦 [2/10] Backing up H100 vectors directory..."
aws s3 sync h100_vectors/ s3://$BUCKET/$PREFIX/h100_vectors/ \
    --no-progress &
PID2=$!

# 3. Backup H100 flat indexes (148 GB)
echo "📦 [3/10] Backing up H100 flat indexes directory..."
aws s3 sync h100_flat_indexes/ s3://$BUCKET/$PREFIX/h100_flat_indexes/ \
    --no-progress &
PID3=$!

# 4. Backup all Python scripts
echo "📦 [4/10] Backing up all Python scripts..."
aws s3 sync . s3://$BUCKET/$PREFIX/scripts/ \
    --exclude "*" \
    --include "*.py" \
    --no-progress &
PID4=$!

# 5. Backup all logs
echo "📦 [5/10] Backing up all log files..."
aws s3 sync . s3://$BUCKET/$PREFIX/logs/ \
    --exclude "*" \
    --include "*.log" \
    --no-progress &
PID5=$!

# 6. Backup all markdown documentation
echo "📦 [6/10] Backing up all documentation..."
aws s3 sync . s3://$BUCKET/$PREFIX/docs/ \
    --exclude "*" \
    --include "*.md" \
    --include "*.txt" \
    --exclude ".venv/*" \
    --no-progress &
PID6=$!

# 7. Backup shell scripts
echo "📦 [7/10] Backing up shell scripts..."
aws s3 sync . s3://$BUCKET/$PREFIX/scripts/ \
    --exclude "*" \
    --include "*.sh" \
    --no-progress &
PID7=$!

# 8. Backup FAISS indexes
echo "📦 [8/10] Backing up FAISS indexes..."
aws s3 sync . s3://$BUCKET/$PREFIX/faiss_indexes/ \
    --exclude "*" \
    --include "*.faiss" \
    --include "*.idx" \
    --include "*.index" \
    --no-progress &
PID8=$!

# 9. Backup numpy arrays and position files
echo "📦 [9/10] Backing up numpy arrays and position files..."
aws s3 sync . s3://$BUCKET/$PREFIX/data/ \
    --exclude "*" \
    --include "*.npy" \
    --include "*positions.txt" \
    --no-progress &
PID9=$!

# 10. Backup reference genome
echo "📦 [10/10] Backing up reference genome..."
aws s3 cp GRCh38.fa s3://$BUCKET/$PREFIX/reference/GRCh38.fa --no-progress &
PID10=$!

echo ""
echo "Waiting for all uploads to complete..."
echo "(This may take several hours for ~300GB+ of data)"
echo ""

# Wait for all background jobs
wait $PID1 && echo "  ✓ Models uploaded"
wait $PID2 && echo "  ✓ H100 vectors uploaded"
wait $PID3 && echo "  ✓ H100 indexes uploaded"
wait $PID4 && echo "  ✓ Python scripts uploaded"
wait $PID5 && echo "  ✓ Logs uploaded"
wait $PID6 && echo "  ✓ Documentation uploaded"
wait $PID7 && echo "  ✓ Shell scripts uploaded"
wait $PID8 && echo "  ✓ FAISS indexes uploaded"
wait $PID9 && echo "  ✓ Numpy arrays uploaded"
wait $PID10 && echo "  ✓ Reference genome uploaded"

echo ""
echo "================================================================================"
echo "✅ COMPREHENSIVE BACKUP COMPLETE!"
echo "================================================================================"
echo ""
echo "Backup location: s3://$BUCKET/$PREFIX/"
echo ""
echo "Folder structure:"
echo "  - models/           All .pt checkpoint files (415 files)"
echo "  - h100_vectors/     All H100 encoded vectors (~148 GB)"
echo "  - h100_flat_indexes/ All FAISS flat indexes (~148 GB)"
echo "  - scripts/          All Python and shell scripts"
echo "  - logs/             All training and validation logs"
echo "  - docs/             All markdown documentation"
echo "  - faiss_indexes/    All FAISS index files"
echo "  - data/             All numpy arrays and position files"
echo "  - reference/        Reference genome (GRCh38.fa)"
echo ""
echo "To list all files:"
echo "  export AWS_ACCESS_KEY_ID='NAKIADXF5GT3EIV7XGCC'"
echo "  export AWS_SECRET_ACCESS_KEY='b2hCHc4rpq+/HWqcF2tB61ATRRRlA9Uf34qEUICt'"
echo "  aws s3 ls s3://$BUCKET/$PREFIX/ --recursive --human-readable --summarize"
echo ""
echo "To download everything:"
echo "  aws s3 sync s3://$BUCKET/$PREFIX/ ./genocache_restored/"
echo ""
echo "Backup timestamp: $TIMESTAMP"
