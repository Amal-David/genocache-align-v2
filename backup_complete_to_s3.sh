#!/bin/bash
# Complete backup of entire genocache directory to S3
# Using genocache-v6 AWS profile

BUCKET="genocache"
PREFIX="genocache_v6.1"
PROFILE="genocache-v6"

echo "================================================================================"
echo "COMPLETE GENOCACHE BACKUP TO S3"
echo "================================================================================"
echo "Bucket: s3://$BUCKET/$PREFIX/"
echo "Profile: $PROFILE"
echo ""

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
echo "Backup timestamp: $TIMESTAMP"
echo ""

# Create a log file for this backup
LOG_FILE="backup_complete_${TIMESTAMP}.log"
echo "Logging to: $LOG_FILE"
echo ""

# Function to log and execute
log_and_run() {
    echo "$1" | tee -a "$LOG_FILE"
    eval "$1" 2>&1 | tee -a "$LOG_FILE"
    return ${PIPESTATUS[0]}
}

# Start the complete sync - this will preserve the entire directory structure
echo "🚀 Starting complete backup of /home/nebius/genocache/ to S3..."
echo "This will preserve all folder structures and upload everything as-is."
echo ""

# Use aws s3 sync to upload everything, excluding only .git directories and temporary files
log_and_run "aws s3 sync /home/nebius/genocache/ s3://$BUCKET/$PREFIX/ \
    --profile $PROFILE \
    --exclude '.git/*' \
    --exclude '*.tmp' \
    --exclude '*.temp' \
    --exclude '.DS_Store'"

if [ $? -eq 0 ]; then
    echo ""
    echo "================================================================================"
    echo "✅ COMPLETE BACKUP SUCCESSFUL!"
    echo "================================================================================"
    echo ""
    echo "Backup location: s3://$BUCKET/$PREFIX/"
    echo "Backup timestamp: $TIMESTAMP"
    echo "Log file: $LOG_FILE"
    echo ""
    echo "To list all files:"
    echo "  aws s3 ls s3://$BUCKET/$PREFIX/ --recursive --human-readable --summarize --profile $PROFILE"
    echo ""
    echo "To download everything:"
    echo "  aws s3 sync s3://$BUCKET/$PREFIX/ ./genocache_restored/ --profile $PROFILE"
    echo ""
else
    echo ""
    echo "❌ BACKUP FAILED!"
    echo "Check the log file: $LOG_FILE"
    echo ""
    exit 1
fi
