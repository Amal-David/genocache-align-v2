#!/bin/bash
#
# GenoCache V4 - Complete Benchmark Pipeline
#
# This script:
# 1. Builds production index (all primary chromosomes)
# 2. Validates the index
# 3. Runs head-to-head benchmark vs minimap2
# 4. Generates summary report
#

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="/home/nebius/genocache/.venv"
LOG_DIR="/home/nebius/genocache/genocache-v4/benchmark/logs"

mkdir -p "$LOG_DIR"

echo "================================================================================"
echo "GenoCache V4 - Complete Benchmark Pipeline"
echo "================================================================================"
echo ""
echo "This will take approximately 1-2 hours:"
echo "  - Index building: ~45-60 minutes (encoding all primary chromosomes)"
echo "  - Benchmark: ~5-10 minutes (testing + minimap2 comparison)"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

# Activate virtual environment
source "$VENV/bin/activate"

# Step 1: Build production index
echo ""
echo "================================================================================"
echo "STEP 1: Building Production Index"
echo "================================================================================"
echo ""

python3 "$SCRIPT_DIR/build_production_index.py" 2>&1 | tee "$LOG_DIR/build_index.log"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "❌ Index building failed! Check logs: $LOG_DIR/build_index.log"
    exit 1
fi

# Step 2: Run benchmark
echo ""
echo "================================================================================"
echo "STEP 2: Running Head-to-Head Benchmark vs minimap2"
echo "================================================================================"
echo ""

python3 "$SCRIPT_DIR/benchmark_vs_minimap2.py" 2>&1 | tee "$LOG_DIR/benchmark.log"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "❌ Benchmark failed! Check logs: $LOG_DIR/benchmark.log"
    exit 1
fi

# Step 3: Summary
echo ""
echo "================================================================================"
echo "BENCHMARK COMPLETE!"
echo "================================================================================"
echo ""
echo "Results:"
echo "  - Index: /home/nebius/genocache/genocache-v4/indexes/"
echo "  - Benchmark: /home/nebius/genocache/genocache-v4/benchmark/results/"
echo "  - Logs: $LOG_DIR/"
echo ""
echo "View results:"
echo "  cat /home/nebius/genocache/genocache-v4/benchmark/results/benchmark_results.json | python3 -m json.tool"
echo ""
