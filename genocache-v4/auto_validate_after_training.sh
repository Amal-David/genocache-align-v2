#!/bin/bash
# Auto-validation script - runs after training completes

LOG_FILE="/home/nebius/genocache/genocache-v4/logs/train_fullgenome.out"
VALIDATION_SCRIPT="/home/nebius/genocache/genocache-v4/validation/validate_fullgenome_pipeline.py"
VENV="/home/nebius/genocache/.venv"

echo "========================================================================"
echo "Auto-Validation Monitor - Waiting for training completion..."
echo "========================================================================"
echo ""

# Wait for training to complete
while true; do
    if grep -q "Training Complete!" "$LOG_FILE" 2>/dev/null; then
        echo "✅ Training completed! Starting validation..."
        sleep 5
        break
    fi
    
    # Check if process is still running
    if ! ps aux | grep -q "[t]rain_fullgenome.py"; then
        echo "⚠️  Training process not found. Checking if completed..."
        if grep -q "Training Complete!" "$LOG_FILE" 2>/dev/null; then
            echo "✅ Training completed! Starting validation..."
            break
        else
            echo "❌ Training process died unexpectedly!"
            exit 1
        fi
    fi
    
    sleep 60  # Check every minute
done

# Run validation
echo ""
echo "========================================================================"
echo "Running Full Genome Pipeline Validation..."
echo "========================================================================"
echo ""

cd /home/nebius/genocache/genocache-v4/validation
source "$VENV/bin/activate"
python3 "$VALIDATION_SCRIPT" 2>&1 | tee /home/nebius/genocache/genocache-v4/logs/auto_validation.log

echo ""
echo "========================================================================"
echo "Validation Complete!"
echo "========================================================================"
echo ""
echo "Results saved to:"
echo "  - /home/nebius/genocache/genocache-v4/logs/auto_validation.log"
echo "  - /home/nebius/genocache/genocache-v4/validation/results/"
