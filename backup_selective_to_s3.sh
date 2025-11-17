#!/bin/bash
# Selective backup of H100 folders and genocache-v4.1-production to S3

BUCKET="genocache"
PREFIX="genocache_v6.1"
PROFILE="genocache-v6"

echo "================================================================================"
echo "SELECTIVE BACKUP TO S3"
echo "================================================================================"
echo "Bucket: s3://$BUCKET/$PREFIX/"
echo "Profile: $PROFILE"
echo ""

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
echo "Backup timestamp: $TIMESTAMP"
echo ""

# Create a log file for this backup
LOG_FILE="backup_selective_${TIMESTAMP}.log"
echo "Logging to: $LOG_FILE"
echo ""

# Function to log and execute
log_and_run() {
    echo "$1" | tee -a "$LOG_FILE"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

echo "📦 Backing up selected folders:"
echo "  - h100_vectors/"
echo "  - h100_flat_indexes/"
echo "  - genocache-v4.1-production/"
echo ""

# Backup h100_vectors
if [ -d "h100_vectors" ]; then
    echo "🚀 [1/3] Backing up h100_vectors..."
    log_and_run "aws s3 sync h100_vectors/ s3://$BUCKET/$PREFIX/h100_vectors/ --profile $PROFILE"
    if [ $? -eq 0 ]; then
        echo "  ✅ h100_vectors backup completed"
    else
        echo "  ❌ h100_vectors backup failed"
        exit 1
    fi
else
    echo "  ⚠️  h100_vectors directory not found"
fi

echo ""

# Backup h100_flat_indexes
if [ -d "h100_flat_indexes" ]; then
    echo "🚀 [2/3] Backing up h100_flat_indexes..."
    log_and_run "aws s3 sync h100_flat_indexes/ s3://$BUCKET/$PREFIX/h100_flat_indexes/ --profile $PROFILE"
    if [ $? -eq 0 ]; then
        echo "  ✅ h100_flat_indexes backup completed"
    else
        echo "  ❌ h100_flat_indexes backup failed"
        exit 1
    fi
else
    echo "  ⚠️  h100_flat_indexes directory not found"
fi

echo ""

# Backup genocache-v4.1-production
if [ -d "genocache-v4.1-production" ]; then
    echo "🚀 [3/3] Backing up genocache-v4.1-production..."
    log_and_run "aws s3 sync genocache-v4.1-production/ s3://$BUCKET/$PREFIX/genocache-v4.1-production/ --profile $PROFILE"
    if [ $? -eq 0 ]; then
        echo "  ✅ genocache-v4.1-production backup completed"
    else
        echo "  ❌ genocache-v4.1-production backup failed"
        exit 1
    fi
else
    echo "  ⚠️  genocache-v4.1-production directory not found"
fi

echo ""
echo "================================================================================"
echo "✅ SELECTIVE BACKUP COMPLETE!"
echo "================================================================================"
echo ""
echo "Backup location: s3://$BUCKET/$PREFIX/"
echo "Backup timestamp: $TIMESTAMP"
echo "Log file: $LOG_FILE"
echo ""
echo "Folders backed up:"
echo "  - s3://$BUCKET/$PREFIX/h100_vectors/"
echo "  - s3://$BUCKET/$PREFIX/h100_flat_indexes/"
echo "  - s3://$BUCKET/$PREFIX/genocache-v4.1-production/"
echo ""
echo "To list all files:"
echo "  aws s3 ls s3://$BUCKET/$PREFIX/ --recursive --human-readable --summarize --profile $PROFILE"
echo ""
