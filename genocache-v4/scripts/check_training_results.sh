#!/bin/bash
# Quick script to check training results

echo "======================================================================"
echo "GenoCache V4 - Training Results Checker"
echo "======================================================================"
echo ""

# Check if training is still running
if ps aux | grep "train_simple.py" | grep -v grep > /dev/null; then
    echo "⏳ Training is STILL RUNNING"
    echo ""
    echo "Current progress:"
    tail -20 logs/training_initial.log
    echo ""
    echo "To monitor live:"
    echo "  tail -f logs/training_initial.log"
else
    echo "✅ Training COMPLETED"
    echo ""
    echo "Final results:"
    tail -50 logs/training_initial.log | grep -A 10 "Epoch 5/5"
    echo ""
    
    # Check if checkpoint exists
    if [ -f "models/checkpoints/best_model.pt" ]; then
        echo "✅ Model checkpoint saved:"
        ls -lh models/checkpoints/best_model.pt
        echo ""
        echo "📊 Key metrics to look for:"
        echo "  - Positive similarity: Should be >0.8"
        echo "  - Negative similarity: Should be <0.3"
        echo "  - Separation: Should be >0.5"
        echo ""
        echo "Next steps:"
        echo "  1. If separation >0.5: Scale to 10M examples"
        echo "  2. If separation >0.7: Ready for curriculum training"
        echo "  3. Check WHAT_TO_DO_NEXT.md for detailed instructions"
    else
        echo "⚠️ No checkpoint found - training may have failed"
    fi
fi

echo ""
echo "======================================================================"
