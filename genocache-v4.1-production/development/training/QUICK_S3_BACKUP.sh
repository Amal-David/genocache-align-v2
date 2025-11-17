#!/bin/bash
# Quick S3 Backup Script - Replace YOUR-BUCKET-NAME with your actual bucket

BUCKET_NAME="YOUR-BUCKET-NAME"  # <<< CHANGE THIS

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                    S3 Backup - Chr22 NAL Training                            ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

if [ "$BUCKET_NAME" = "YOUR-BUCKET-NAME" ]; then
    echo "❌ ERROR: Please edit this script and set BUCKET_NAME to your actual S3 bucket!"
    echo ""
    echo "Edit line 4 of this script:"
    echo '  BUCKET_NAME="your-actual-bucket-name"'
    echo ""
    exit 1
fi

echo "Bucket: s3://$BUCKET_NAME/nal_chr22_backup/"
echo ""
echo "Step 1: Backing up code and documentation (~104 KB)..."
aws s3 sync nal_chr22_PRODUCTION_BACKUP/ \
  s3://$BUCKET_NAME/nal_chr22_backup/ \
  --exclude "models/*" \
  --exclude "indexes/*" \
  --exclude "__pycache__/*" \
  --exclude "*.pyc"

echo ""
echo "Step 2: Verifying backup..."
aws s3 ls s3://$BUCKET_NAME/nal_chr22_backup/ --recursive --human-readable

echo ""
echo "✅ Backup complete!"
echo ""
echo "To restore:"
echo "  aws s3 sync s3://$BUCKET_NAME/nal_chr22_backup/ ./restored/"
echo ""
echo "Optional: Backup models separately (~5.5 MB)"
echo "  aws s3 cp nal_single_chr/models/chr22_nal_512bp_final.pt s3://$BUCKET_NAME/nal_chr22_backup/models/"
