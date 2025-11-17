#!/bin/bash
# Quick status check - run this anytime after reconnecting to VM

echo "========================================"
echo "GenoCache V3 Training Status"
echo "========================================"
echo ""

# Check if PID file exists
if [ -f training.pid ]; then
    PID=$(cat training.pid)
    
    # Check if process is still running
    if ps -p $PID > /dev/null 2>&1; then
        echo "✅ Training is RUNNING"
        echo "   PID: $PID"
        
        # Get elapsed time
        ELAPSED=$(ps -p $PID -o etime= | tr -d ' ')
        echo "   Running for: $ELAPSED"
        
        # GPU status
        echo ""
        echo "GPU Status:"
        nvidia-smi --query-gpu=utilization.gpu,memory.used,temperature.gpu --format=csv,noheader,nounits | \
        awk -F', ' '{printf "   Utilization: %s%%\n   Memory: %s MB\n   Temp: %s°C\n", $1, $2, $3}'
    else
        echo "❌ Training process STOPPED (PID $PID not found)"
        echo "   Check training.log for errors"
    fi
else
    echo "⚠️  No PID file found"
    echo "   Searching for training process..."
    
    RUNNING_PID=$(pgrep -f "train_genocache.py" | head -1)
    if [ -n "$RUNNING_PID" ]; then
        echo "   Found: PID $RUNNING_PID"
        echo $RUNNING_PID > training.pid
    else
        echo "   No training process found"
    fi
fi

echo ""
echo "========================================"
echo "Recent Log Output:"
echo "========================================"
tail -20 training.log

echo ""
echo "========================================"
echo "Models Saved:"
echo "========================================"
ls -lh models/*.pt 2>/dev/null || echo "   No models saved yet"

echo ""
echo "Commands:"
echo "  Full log: tail -f training.log"
echo "  Stop: kill \$(cat training.pid)"
echo "  Monitor: ./monitor_training.sh"
echo ""
