#!/bin/bash
# Backup important model files to S3

# S3 Configuration
export AWS_ACCESS_KEY_ID="NAKIADXF5GT3EIV7XGCC"
export AWS_SECRET_ACCESS_KEY="b2hCHc4rpq+/HWqcF2tB61ATRRRlA9Uf34qEUICt"
export AWS_DEFAULT_REGION="us-east-1"
BUCKET="hackathon-team-fabric3-5"
PREFIX="genocache_models"

echo "================================================================================"
echo "BACKING UP GENOCACHE MODELS TO S3"
echo "================================================================================"
echo "Bucket: s3://$BUCKET/$PREFIX/"
echo ""

# Create timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
echo "Backup timestamp: $TIMESTAMP"
echo ""

# Backup best models (most important)
echo "📦 Backing up BEST models..."
aws s3 cp cnn_hyenadna_large_best.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp experimental_model_1_cnn_transformer.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp experimental_model_2_improved_cnn_rnn.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp h100_model_1_cnn.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp h100_model_2_deep_cnn.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp h100_model_3_cnn_rnn.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp h100_model_4_transformer.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress
aws s3 cp h100_model_5_dilated_cnn.pt s3://$BUCKET/$PREFIX/best_models/ --no-progress

echo ""
echo "📦 Backing up CNN-HyenaDNA epoch checkpoints..."
aws s3 sync . s3://$BUCKET/$PREFIX/cnn_hyenadna_epochs/ \
    --exclude "*" \
    --include "cnn_hyenadna_large_epoch*.pt" \
    --no-progress

echo ""
echo "📦 Backing up training scripts..."
aws s3 cp train_cnn_hyenadna.py s3://$BUCKET/$PREFIX/scripts/ --no-progress
aws s3 cp train_ensemble_models.py s3://$BUCKET/$PREFIX/scripts/ --no-progress
aws s3 cp train_experimental_models.py s3://$BUCKET/$PREFIX/scripts/ --no-progress
aws s3 cp encode_cnn_hyenadna.py s3://$BUCKET/$PREFIX/scripts/ --no-progress

echo ""
echo "📦 Backing up training logs..."
aws s3 cp cnn_hyenadna_training_gradaccum.log s3://$BUCKET/$PREFIX/logs/ --no-progress
aws s3 cp experimental_training.log s3://$BUCKET/$PREFIX/logs/ --no-progress
aws s3 cp h100_training.log s3://$BUCKET/$PREFIX/logs/ --no-progress

echo ""
echo "📦 Backing up analysis documents..."
aws s3 sync . s3://$BUCKET/$PREFIX/docs/ \
    --exclude "*" \
    --include "*.md" \
    --no-progress

echo ""
echo "================================================================================"
echo "✅ BACKUP COMPLETE!"
echo "================================================================================"
echo ""
echo "Uploaded to: s3://$BUCKET/$PREFIX/"
echo ""
echo "To list files:"
echo "  aws s3 ls s3://$BUCKET/$PREFIX/ --recursive"
echo ""
echo "To download:"
echo "  aws s3 sync s3://$BUCKET/$PREFIX/ ./restored_models/"
