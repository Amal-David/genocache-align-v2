#!/bin/bash
# Automated pipeline: wait for data generation → train → validate

set -e

echo "================================================================================"
echo "GenoCache V4 - Automated Pipeline"
echo "================================================================================"
echo ""

cd /home/nebius/genocache/genocache-v4

# Step 1: Wait for data generation to complete
echo "Step 1: Waiting for data generation to complete..."
while ps aux | grep "generate_training_data_10M.py" | grep -v grep > /dev/null; do
    sleep 30
    tail -2 logs/generate_10M.log | grep "Progress:" || true
done

echo "✅ Data generation complete!"
echo ""

# Check if data file exists
if [ ! -f "data/training_10M.h5" ]; then
    echo "❌ ERROR: training_10M.h5 not found!"
    exit 1
fi

FILE_SIZE=$(du -h data/training_10M.h5 | cut -f1)
echo "✅ Data file created: $FILE_SIZE"
echo ""

# Step 2: Start training
echo "================================================================================"
echo "Step 2: Starting training on 10M examples..."
echo "================================================================================"
echo ""

.venv/bin/python training/train_10M.py > logs/training_10M.log 2>&1

echo "✅ Training complete!"
echo ""

# Step 3: Run validation tests
echo "================================================================================"
echo "Step 3: Running validation tests..."
echo "================================================================================"
echo ""

.venv/bin/python scripts/test_model_10M.py > logs/validation_10M.log 2>&1

echo "✅ Validation complete!"
echo ""

# Step 4: Summary
echo "================================================================================"
echo "Pipeline Complete - Results Summary"
echo "================================================================================"
echo ""

echo "Training results:"
tail -30 logs/training_10M.log | grep -E "Best separation|Epoch 10/10|Results" | tail -10

echo ""
echo "Validation results:"
tail -20 logs/validation_10M.log

echo ""
echo "================================================================================"
echo "✅ Full pipeline complete!"
echo "================================================================================"
