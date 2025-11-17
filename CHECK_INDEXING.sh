#!/bin/bash
# Monitor NAL indexing progress

cd /home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                      NAL INDEXING STATUS                                     ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check if process is running
PID=$(cat indexing.pid 2>/dev/null)
if [ -n "$PID" ] && ps -p $PID > /dev/null 2>&1; then
    echo "✅ Indexing is RUNNING"
    echo "   PID: $PID"
    ps -p $PID -o pid,etime,%cpu,%mem,cmd --no-headers
    echo ""
else
    echo "❌ Indexing is NOT running"
    echo ""
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " LATEST OUTPUT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

tail -30 indexing_run.log 2>/dev/null || tail -30 logs/indexing_*.log 2>/dev/null | tail -30

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check if index file exists
if [ -f "indexes/genocache_nal_stride32.index" ]; then
    echo "✅ Index file exists: indexes/genocache_nal_stride32.index"
    ls -lh indexes/genocache_nal_stride32.index | awk '{print "   Size:", $5, " Modified:", $6, $7, $8}'
    
    if [ -f "indexes/genocache_nal_stride32.positions.npz" ]; then
        echo "✅ Positions file exists: indexes/genocache_nal_stride32.positions.npz"
        ls -lh indexes/genocache_nal_stride32.positions.npz | awk '{print "   Size:", $5, " Modified:", $6, $7, $8}'
    fi
else
    echo "⏳ Index file not yet created"
fi

echo ""
echo "To monitor live: watch -n 10 /home/nebius/genocache/CHECK_INDEXING.sh"
