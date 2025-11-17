#!/bin/bash
#
# Test GenoCache V4 Complete Pipeline
# Tests adaptive seeding + WFA alignment
#

set -e

VENV="/home/nebius/genocache/.venv"
SCRIPT_DIR="/home/nebius/genocache/genocache-v4"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║      GenoCache V4 - Test Complete Pipeline                         ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

# Activate venv
source "$VENV/bin/activate"

# Check files exist
echo "Checking required files..."
echo ""

MODEL="$SCRIPT_DIR/models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt"
INDEX="$SCRIPT_DIR/indexes/genocache_v4_production.index"
METADATA="$SCRIPT_DIR/indexes/genocache_v4_production.metadata.pkl"
GENOME="/home/nebius/genocache/GRCh38.fa"
READS="/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa"

if [ ! -f "$MODEL" ]; then
    echo "❌ Model not found: $MODEL"
    exit 1
fi
echo "✅ Model found"

if [ ! -f "$INDEX" ]; then
    echo "❌ Index not found: $INDEX"
    exit 1
fi
echo "✅ Index found"

if [ ! -f "$METADATA" ]; then
    echo "❌ Metadata not found: $METADATA"
    exit 1
fi
echo "✅ Metadata found"

if [ ! -f "$GENOME" ]; then
    echo "❌ Genome not found: $GENOME"
    exit 1
fi
echo "✅ Genome found"

if [ ! -f "$READS" ]; then
    echo "❌ Reads not found: $READS"
    exit 1
fi
echo "✅ Reads found ($(grep -c '^>' $READS) reads)"
echo ""

# Test on small subset first (10 reads)
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 1: Small test (10 reads)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python3 "$SCRIPT_DIR/complete_pipeline.py" \
    --model "$MODEL" \
    --index "$INDEX" \
    --metadata "$METADATA" \
    --genome "$GENOME" \
    --reads "$READS" \
    --output "$SCRIPT_DIR/test_output_10.sam" \
    --max-reads 10 \
    --device cuda

echo ""
echo "✅ Small test complete!"
echo ""
echo "Output: $SCRIPT_DIR/test_output_10.sam"
echo ""

# Check output
if [ -f "$SCRIPT_DIR/test_output_10.sam" ]; then
    LINES=$(wc -l < "$SCRIPT_DIR/test_output_10.sam")
    echo "SAM file has $LINES lines"
    echo ""
    echo "First few lines:"
    head -20 "$SCRIPT_DIR/test_output_10.sam"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 2: Full test (100 reads)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python3 "$SCRIPT_DIR/complete_pipeline.py" \
    --model "$MODEL" \
    --index "$INDEX" \
    --metadata "$METADATA" \
    --genome "$GENOME" \
    --reads "$READS" \
    --output "$SCRIPT_DIR/test_output_100.sam" \
    --max-reads 100 \
    --device cuda

echo ""
echo "✅ Full test complete!"
echo ""
echo "Output: $SCRIPT_DIR/test_output_100.sam"
echo ""

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    TESTS COMPLETE! ✅                               ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Next steps:"
echo "  1. Test on real GIAB data"
echo "  2. Compare with minimap2 on same data"
echo "  3. Validate accuracy and speed"
echo ""
