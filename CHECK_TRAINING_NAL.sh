#!/bin/bash
# Quick script to monitor NAL training progress

cd /home/nebius/genocache/genocache-v4.1-production

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                      NAL TRAINING STATUS                                     ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check if process is running
PID=$(cat ../nal_training.pid 2>/dev/null)
if [ -n "$PID" ] && ps -p $PID > /dev/null 2>&1; then
    echo "✅ Training is RUNNING"
    echo "   PID: $PID"
    ps -p $PID -o pid,etime,%cpu,%mem,cmd --no-headers
    echo ""
else
    echo "❌ Training is NOT running"
    echo ""
fi

# Show latest log entries
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " LATEST TRAINING METRICS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

LOG_FILE=$(ls -t models/training_nal_*.log 2>/dev/null | head -1)
if [ -f "$LOG_FILE" ]; then
    echo "Log: $LOG_FILE"
    echo ""
    echo "Epoch | Train Loss | Val Loss | Pos Sim | Neg Sim | Separation"
    echo "------|------------|----------|---------|---------|------------"
    tail -10 "$LOG_FILE" | awk -F',' '{if (NR>1) printf "%-6s| %-11.4f| %-9.4f| %-8.3f| %-8.3f| %.3f\n", $1, $2, $6, $3, $4, $5}'
else
    echo "No log file found yet"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check model file
if [ -f "models/genocache_nal.pt" ]; then
    echo "✅ Model checkpoint exists: models/genocache_nal.pt"
    ls -lh models/genocache_nal.pt | awk '{print "   Size:", $5, " Modified:", $6, $7, $8}'
else
    echo "⏳ Model checkpoint not yet created"
fi

echo ""
echo "To monitor live: watch -n 5 ./CHECK_TRAINING_NAL.sh"
