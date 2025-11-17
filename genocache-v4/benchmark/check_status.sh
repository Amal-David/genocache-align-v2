#!/bin/bash
# Quick status checker for benchmark pipeline

LOG_FILE="/home/nebius/genocache/genocache-v4/benchmark/logs/full_benchmark.out"
BUILD_LOG="/home/nebius/genocache/genocache-v4/benchmark/logs/build_index.log"
BENCH_LOG="/home/nebius/genocache/genocache-v4/benchmark/logs/benchmark.log"
PID_FILE="/home/nebius/genocache/genocache-v4/benchmark/logs/benchmark.pid"

echo "========================================================================"
echo "GenoCache V4 Benchmark - Status Check"
echo "========================================================================"
echo ""

# Check if process is running
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if ps -p $PID > /dev/null 2>&1; then
        echo "✅ Benchmark is RUNNING (PID: $PID)"
        
        # Determine current phase
        if grep -q "STEP 2:" "$LOG_FILE" 2>/dev/null; then
            echo "📍 Phase: Running benchmark vs minimap2"
        elif grep -q "STEP 1:" "$LOG_FILE" 2>/dev/null; then
            echo "📍 Phase: Building production index"
            
            # Show encoding progress
            if [ -f "$BUILD_LOG" ]; then
                echo ""
                echo "Encoding progress:"
                tail -5 "$BUILD_LOG" | grep -E "NC_|chr" | tail -3
            fi
        fi
    else
        echo "⚠️  Process not running (PID $PID no longer exists)"
        
        # Check if completed
        if grep -q "BENCHMARK COMPLETE" "$LOG_FILE" 2>/dev/null; then
            echo "✅ Benchmark COMPLETED!"
        else
            echo "❌ Benchmark may have failed. Check logs."
        fi
    fi
else
    echo "❌ No benchmark running (no PID file found)"
fi

echo ""
echo "========================================================================"
echo "Quick Logs:"
echo "========================================================================"
tail -20 "$LOG_FILE" 2>/dev/null || echo "No logs yet"

echo ""
echo "========================================================================"
echo "Monitor live:"
echo "  tail -f $LOG_FILE"
echo ""
echo "GPU status:"
nvidia-smi --query-gpu=utilization.gpu,memory.used,temperature.gpu --format=csv,noheader
