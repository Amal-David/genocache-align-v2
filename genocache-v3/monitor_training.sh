#!/bin/bash
# Monitor GenoCache V3 training progress

echo "================================"
echo "GenoCache V3 Training Monitor"
echo "================================"
echo ""

# Check if training is running
PID=$(pgrep -f "train_genocache.py")
if [ -z "$PID" ]; then
    echo "❌ Training process not found"
else
    echo "✅ Training running (PID: $PID)"
    
    # GPU utilization
    echo ""
    echo "GPU Status:"
    nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits | \
    awk -F', ' '{printf "  GPU %s: %s\n  Utilization: %s%%\n  Memory: %s / %s MB\n", $1, $2, $3, $4, $5}'
fi

echo ""
echo "================================"
echo "Training Log (last 30 lines):"
echo "================================"
tail -30 training.log

echo ""
echo "================================"
echo "Commands:"
echo "================================"
echo "  Watch live: tail -f training.log"
echo "  Stop training: kill $PID"
echo "  Check models: ls -lh models/"
echo ""
