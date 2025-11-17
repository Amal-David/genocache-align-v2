#!/bin/bash
# Test on 1000 reads to validate scalability and avoid overfitting

set -e

cd "$(dirname "$0")"
source /home/nebius/genocache/.venv/bin/activate

BASE_DIR="/home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned"
MODEL="${BASE_DIR}/models/genocache_nal.pt"
INDEX="${BASE_DIR}/indexes/genocache_nal_stride32.index"
POSITIONS="${BASE_DIR}/indexes/genocache_nal_stride32.positions.npz"
REFERENCE="/home/nebius/genocache/GRCh38.fa"

# Original full GIAB dataset
GIAB_FULL="/home/nebius/genocache/giab_hg002_chr22_ont.fastq.gz"

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                    1000 READ VALIDATION TEST                                 ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Purpose: Validate that optimization didn't overfit to 100-read test set"
echo ""

# Extract 1000 reads (skip first 100 we used for testing)
if [ ! -f "data/giab_1000reads.fastq" ]; then
    echo "Extracting 1000 reads from GIAB dataset..."
    mkdir -p data
    
    if [ -f "$GIAB_FULL" ]; then
        # Skip first 100 reads, take next 1000
        zcat "$GIAB_FULL" | awk 'BEGIN {skip=100; count=0; take=1000} 
            NR % 4 == 1 {
                if (skip > 0) {
                    skip--;
                    next;
                }
                if (count >= take) {
                    exit;
                }
                count++;
            }
            skip == 0 && count <= take' > data/giab_1000reads.fastq
        
        actual=$(grep -c "^@" data/giab_1000reads.fastq || echo 0)
        echo "  ✅ Extracted $actual reads"
    else
        echo "  ⚠️  Full GIAB file not found: $GIAB_FULL"
        echo "  Using alternative: taking 1000 reads from any available dataset"
        
        # Try to find any FASTQ file
        FASTQ_FILE=$(find /home/nebius/genocache -name "*.fastq" -o -name "*.fastq.gz" | head -1)
        if [ -n "$FASTQ_FILE" ]; then
            if [[ "$FASTQ_FILE" == *.gz ]]; then
                zcat "$FASTQ_FILE" | head -4000 > data/giab_1000reads.fastq
            else
                head -4000 "$FASTQ_FILE" > data/giab_1000reads.fastq
            fi
            echo "  ✅ Extracted from: $FASTQ_FILE"
        else
            echo "  ❌ No FASTQ files found!"
            exit 1
        fi
    fi
else
    echo "Using existing data/giab_1000reads.fastq"
fi

READS_1000="data/giab_1000reads.fastq"
actual_count=$(grep -c "^@" "$READS_1000")
echo "  Total reads in dataset: $actual_count"
echo ""

# Test 1: minimap2 baseline
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 1: minimap2 Baseline"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if command -v minimap2 >/dev/null 2>&1; then
    echo "Running minimap2..."
    time minimap2 -ax map-ont "$REFERENCE" "$READS_1000" > results/minimap2_1000reads.sam 2>&1
    
    mm2_total=$(grep -v "^@" results/minimap2_1000reads.sam | awk '!and($2, 256) && !and($2, 2048)' | wc -l)
    mm2_mapped=$(grep -v "^@" results/minimap2_1000reads.sam | awk '!and($2, 256) && !and($2, 2048) && $3 != "*"' | wc -l)
    mm2_rate=$(echo "scale=1; $mm2_mapped * 100 / $mm2_total" | bc)
    echo "  minimap2: $mm2_mapped/$mm2_total = $mm2_rate%"
else
    echo "  ⚠️  minimap2 not found, skipping baseline"
fi

# Test 2: 3-Tier NAL (user's proposed config)
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 2: NAL 3-Tier (6→12→32, K=64/64/48)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "Running NAL 3-tier on 1000 reads..."
time python3 align_adaptive_3tier.py \
    --model "$MODEL" \
    --index "$INDEX" \
    --positions "$POSITIONS" \
    --reference "$REFERENCE" \
    --reads "$READS_1000" \
    --output "results/nal_3tier_1000reads.sam" \
    --device cuda

# Calculate mapping rate
nal_total=$(grep -v "^@" results/nal_3tier_1000reads.sam | wc -l)
nal_mapped=$(grep -v "^@" results/nal_3tier_1000reads.sam | awk '$3 != "*"' | wc -l)
nal_rate=$(echo "scale=1; $nal_mapped * 100 / $nal_total" | bc)

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RESULTS SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "results/minimap2_1000reads.sam" ]; then
    echo "  minimap2:    $mm2_mapped/$mm2_total = $mm2_rate%"
fi
echo "  NAL 3-tier:  $nal_mapped/$nal_total = $nal_rate%"

if [ -f "results/minimap2_1000reads.sam" ]; then
    diff=$(echo "scale=1; $nal_rate - $mm2_rate" | bc)
    echo ""
    echo "  Difference: $diff%"
    
    if (( $(echo "$diff >= 0" | bc -l) )); then
        echo "  ✅ NAL matches or exceeds minimap2!"
    else
        echo "  ⚠️  NAL below minimap2 by $(echo "$diff * -1" | bc)%"
    fi
fi

# False positive detection
if [ -f "results/minimap2_1000reads.sam" ]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "FALSE POSITIVE DETECTION (1000 reads)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    python3 detect_false_positives.py \
        --nal results/nal_3tier_1000reads.sam \
        --minimap2 results/minimap2_1000reads.sam \
        --tolerance 1000 \
        --nal-name "NAL-3tier"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "VALIDATION COMPLETE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Files generated:"
echo "  • data/giab_1000reads.fastq           (test dataset)"
echo "  • results/minimap2_1000reads.sam       (baseline)"
echo "  • results/nal_3tier_1000reads.sam      (NAL results)"
echo ""

if (( $(echo "$nal_rate >= 90" | bc -l) )); then
    echo "✅ VALIDATION PASSED: ≥90% mapping rate on 1000 reads"
    echo "   No significant overfitting detected!"
else
    echo "⚠️  VALIDATION CONCERN: <90% mapping rate"
    echo "   May need parameter adjustment for diverse data"
fi
echo ""

