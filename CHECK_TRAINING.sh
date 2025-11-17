#!/bin/bash
# Quick training status check

echo "════════════════════════════════════════"
echo "  TRAINING STATUS"
echo "════════════════════════════════════════"
echo ""

# Check screen session
if screen -ls | grep -q "full_genome_training"; then
    echo "✓ Screen session: RUNNING"
else
    echo "✗ Screen session: NOT FOUND"
    exit 1
fi

# Check log file
if [ -f full_genome_training.log ]; then
    LATEST_EPOCH=$(grep "Epoch" full_genome_training.log 2>/dev/null | tail -1)
    LATEST_LOSS=$(echo "$LATEST_EPOCH" | grep -oP "Loss = \K[0-9.]+")
    
    echo ""
    echo "Latest progress:"
    echo "$LATEST_EPOCH"
    echo ""
    
    # Count epochs completed
    EPOCHS_DONE=$(grep -c "Epoch.*Loss" full_genome_training.log 2>/dev/null)
    echo "Epochs completed: $EPOCHS_DONE / 20"
    PERCENT=$((EPOCHS_DONE * 100 / 20))
    echo "Progress: $PERCENT%"
    echo ""
    
    # Time estimate
    if [ $EPOCHS_DONE -gt 0 ]; then
        REMAINING=$((20 - EPOCHS_DONE))
        MIN_REMAINING=$((REMAINING * 6))
        MAX_REMAINING=$((REMAINING * 8))
        echo "Estimated time remaining: $MIN_REMAINING-$MAX_REMAINING minutes"
    fi
else
    echo "✗ Log file not found"
fi

echo ""
echo "════════════════════════════════════════"
echo "Commands:"
echo "  screen -r full_genome_training  # Watch live"
echo "  tail -f full_genome_training.log  # Follow log"
echo "════════════════════════════════════════"
