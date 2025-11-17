#!/bin/bash

# NAL Pipeline Testing Script
# Tests complete alignment pipeline on GIAB data

set -e

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                  NAL Pipeline Testing - GIAB Data                             ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
BASE_DIR="/home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned"
MODEL="${BASE_DIR}/models/genocache_nal.pt"
INDEX="${BASE_DIR}/indexes/genocache_nal_stride32.index"
POSITIONS="${BASE_DIR}/indexes/genocache_nal_stride32.positions.npz"
REFERENCE="/home/nebius/genocache/GRCh38.fa"
READS="/home/nebius/genocache/genocache-v4.1-production/validation/data/giab_hg002_100reads.fastq"
OUTPUT="${BASE_DIR}/results/giab_nal_aligned.sam"

echo "Configuration:"
echo "  Model: ${MODEL}"
echo "  Index: ${INDEX}"
echo "  Reference: ${REFERENCE}"
echo "  Reads: ${READS}"
echo "  Output: ${OUTPUT}"
echo ""

# Check files exist
echo "Checking files..."
for file in "${MODEL}" "${INDEX}" "${POSITIONS}" "${REFERENCE}" "${READS}"; do
    if [ ! -f "$file" ]; then
        echo "  ❌ Missing: $file"
        exit 1
    fi
done
echo "  ✅ All files present"
echo ""

# Create results directory
mkdir -p "${BASE_DIR}/results"

# Test individual components first
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1. Testing Seeding Module"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
cd "${BASE_DIR}"
python3 seeding_nal.py "${MODEL}" "${INDEX}" "${POSITIONS}" || {
    echo "❌ Seeding test failed"
    exit 1
}
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2. Testing Chaining Module"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
python3 chaining_nal.py || {
    echo "❌ Chaining test failed"
    exit 1
}
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3. Testing Complete Pipeline on GIAB 100 Reads"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Starting alignment..."
echo ""

python3 align_nal.py \
    --model "${MODEL}" \
    --index "${INDEX}" \
    --positions "${POSITIONS}" \
    --reference "${REFERENCE}" \
    --reads "${READS}" \
    --output "${OUTPUT}" || {
    echo "❌ Pipeline test failed"
    exit 1
}

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4. Analyzing Results"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Count mapped reads
if [ -f "${OUTPUT}" ]; then
    total=$(grep -v "^@" "${OUTPUT}" | wc -l)
    mapped=$(grep -v "^@" "${OUTPUT}" | awk '$3 != "*"' | wc -l)
    unmapped=$(grep -v "^@" "${OUTPUT}" | awk '$3 == "*"' | wc -l)
    
    if [ "$total" -gt 0 ]; then
        mapping_rate=$(echo "scale=1; $mapped * 100 / $total" | bc)
    else
        mapping_rate=0
    fi
    
    echo "Results:"
    echo "  Total reads: ${total}"
    echo "  Mapped: ${mapped} (${mapping_rate}%)"
    echo "  Unmapped: ${unmapped}"
    echo ""
    
    echo "SAM output: ${OUTPUT}"
    echo "  Size: $(du -h "${OUTPUT}" | cut -f1)"
    echo ""
    
    # Show first few alignments
    echo "First 5 mapped reads:"
    grep -v "^@" "${OUTPUT}" | awk '$3 != "*"' | head -5 | \
        awk '{printf "  %s → %s:%s (MAPQ=%s)\n", $1, $3, $4, $5}'
    echo ""
    
    echo "╔══════════════════════════════════════════════════════════════════════════════╗"
    if (( $(echo "$mapping_rate >= 80" | bc -l) )); then
        echo "║                        ✅ SUCCESS! ✅                                        ║"
        echo "╚══════════════════════════════════════════════════════════════════════════════╝"
        echo ""
        echo "Mapping rate: ${mapping_rate}% - Target achieved! 🎉"
    elif (( $(echo "$mapping_rate >= 50" | bc -l) )); then
        echo "║                    ⚠️  PARTIAL SUCCESS ⚠️                                   ║"
        echo "╚══════════════════════════════════════════════════════════════════════════════╝"
        echo ""
        echo "Mapping rate: ${mapping_rate}% - Better than baseline (32%), but below target (85%)"
    else
        echo "║                       ❌ NEEDS IMPROVEMENT ❌                                ║"
        echo "╚══════════════════════════════════════════════════════════════════════════════╝"
        echo ""
        echo "Mapping rate: ${mapping_rate}% - Below expectations"
    fi
else
    echo "❌ Output file not created"
    exit 1
fi

echo ""
echo "Next steps:"
echo "  1. Compare with minimap2: minimap2 -ax map-ont ${REFERENCE} ${READS} > minimap2.sam"
echo "  2. Validate accuracy with ground truth"
echo "  3. Test on larger dataset"
echo ""
