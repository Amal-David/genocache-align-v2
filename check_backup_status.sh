#!/bin/bash
# Monitor backup progress

BUCKET="genocache"
PREFIX="genocache_v6.1"
PROFILE="genocache-v6"

echo "================================================================================"
echo "BACKUP STATUS CHECK"
echo "================================================================================"
echo ""

# Check if backup process is still running
if pgrep -f "backup_complete_to_s3.sh" > /dev/null; then
    echo "✅ Backup process is still running"
    echo ""
    
    # Show latest log entries
    LATEST_LOG=$(ls -t backup_complete_*.log 2>/dev/null | head -1)
    if [ -n "$LATEST_LOG" ]; then
        echo "📋 Latest log entries from $LATEST_LOG:"
        tail -10 "$LATEST_LOG"
        echo ""
    fi
else
    echo "⏹️  Backup process has finished (or not running)"
    echo ""
fi

# Check what's been uploaded to S3
echo "📊 Current S3 bucket contents:"
aws s3 ls s3://$BUCKET/$PREFIX/ --recursive --human-readable --summarize --profile $PROFILE

echo ""
echo "================================================================================"
echo "To check detailed progress, run: tail -f backup_complete_*.log"
echo "================================================================================"
