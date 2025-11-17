#!/bin/bash
# Monitor Chr22 focused training progress

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                    CHR22 FOCUSED TRAINING - PROGRESS                         ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check if training is running
if ps aux | grep -q "[t]rain_chr22_focused.py"; then
    echo "✅ Training is RUNNING"
else
    echo "❌ Training is NOT running"
fi

echo ""
echo "Latest progress:"
echo "────────────────────────────────────────────────────────────────────────────────"
tail -20 logs/training.log | grep "Batch"

echo ""
echo "Current stage:"
LATEST_BATCH=$(tail -100 logs/training.log | grep "Batch" | tail -1 | awk -F'/' '{print $1}' | awk '{print $2}')
if [ ! -z "$LATEST_BATCH" ]; then
    PERCENT=$(echo "scale=1; $LATEST_BATCH * 100 / 8000" | bc)
    echo "  Progress: $LATEST_BATCH/8000 batches ($PERCENT%)"
    
    if [ "$LATEST_BATCH" -lt 1000 ]; then
        echo "  Stage: Warm-up (5% error rate)"
    elif [ "$LATEST_BATCH" -lt 2500 ]; then
        echo "  Stage: Early training (8% error rate)"
    elif [ "$LATEST_BATCH" -lt 5000 ]; then
        echo "  Stage: Mid training (12% error rate)"
    else
        echo "  Stage: Advanced training (15% error rate)"
    fi
fi

echo ""
echo "Saved checkpoints:"
echo "────────────────────────────────────────────────────────────────────────────────"
ls -lh models/*.pt 2>/dev/null | tail -5

echo ""
echo "To view live log: tail -f logs/training.log"
echo ""
