#!/bin/bash
# Compare 2-tier vs 3-tier strategies and detect false positives

set -e

cd "$(dirname "$0")"
source /home/nebius/genocache/.venv/bin/activate

BASE_DIR="/home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned"
MODEL="${BASE_DIR}/models/genocache_nal.pt"
INDEX="${BASE_DIR}/indexes/genocache_nal_stride32.index"
POSITIONS="${BASE_DIR}/indexes/genocache_nal_stride32.positions.npz"
REFERENCE="/home/nebius/genocache/GRCh38.fa"
READS="${BASE_DIR}/../../../validation/data/giab_hg002_100reads.fastq"
MINIMAP2_BASELINE="${BASE_DIR}/results/minimap2_baseline.sam"

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║              TIER STRATEGY COMPARISON & FALSE POSITIVE DETECTION             ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Testing strategies:"
echo "  • 2-Tier: 6→12 seeds (K=64)"
echo "  • 3-Tier: 6→12→32 seeds (K=64, K=64, K=48)"
echo "  • Baseline: minimap2"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Test 1: Three-tier (user's proposal)
echo ""
echo "Test 1: Three-Tier Strategy (6→12→32, K=64/64/48)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

python3 align_adaptive_3tier.py \
    --model "$MODEL" \
    --index "$INDEX" \
    --positions "$POSITIONS" \
    --reference "$REFERENCE" \
    --reads "$READS" \
    --output "results/adaptive_3tier.sam" \
    --device cuda

# Test 2: Two-tier (our recommendation)
echo ""
echo "Test 2: Two-Tier Strategy (6→12, K=64/64)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Modify 3-tier script to skip tier 3 (quick hack)
python3 << 'EOF'
import sys
sys.path.insert(0, '..')
from test.align_adaptive_3tier import *

# Override tiers for 2-tier
aligner = ThreeTierAligner(
    model_path='../models/genocache_nal.pt',
    index_path='../indexes/genocache_nal_stride32.index',
    positions_path='../indexes/genocache_nal_stride32.positions.npz',
    reference_path='/home/nebius/genocache/GRCh38.fa',
    device='cuda'
)

# Override to 2-tier
aligner.tiers = [
    {'name': 'Fast', 'seeds': 6, 'K': 64, 'tolerance': 1000, 'exit_score': 3},
    {'name': 'Standard', 'seeds': 12, 'K': 64, 'tolerance': 1500, 'exit_score': None},
]

# Load reads
reads = load_reads('/home/nebius/genocache/genocache-v4.1-production/validation/data/giab_hg002_100reads.fastq', max_reads=100)
print(f"Loaded {len(reads)} reads")

# Align
print("\nAligning with 2-tier strategy...")
results = aligner.align_reads(reads, verbose=True)

# Write SAM
aligner.write_sam(results, 'results/adaptive_2tier.sam')
print("\n✅ 2-tier alignment complete!")

EOF

# Compare mapping rates
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "MAPPING RATE COMPARISON:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

for strategy in 2tier 3tier; do
    sam="results/adaptive_${strategy}.sam"
    if [ -f "$sam" ]; then
        total=$(grep -v "^@" "$sam" | wc -l)
        mapped=$(grep -v "^@" "$sam" | awk '$3 != "*"' | wc -l)
        rate=$(echo "scale=1; $mapped * 100 / $total" | bc)
        echo "  ${strategy}: $mapped/$total = $rate%"
    fi
done

mm2_total=$(grep -v "^@" "$MINIMAP2_BASELINE" | awk '!and($2, 256) && !and($2, 2048)' | wc -l)
mm2_mapped=$(grep -v "^@" "$MINIMAP2_BASELINE" | awk '!and($2, 256) && !and($2, 2048) && $3 != "*"' | wc -l)
mm2_rate=$(echo "scale=1; $mm2_mapped * 100 / $mm2_total" | bc)
echo "  minimap2: $mm2_mapped/$mm2_total = $mm2_rate%"

# False positive detection
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "FALSE POSITIVE DETECTION:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo ""
echo "Comparing 2-Tier vs minimap2:"
python3 detect_false_positives.py \
    --nal results/adaptive_2tier.sam \
    --minimap2 "$MINIMAP2_BASELINE" \
    --tolerance 1000 \
    --nal-name "NAL-2tier"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Comparing 3-Tier vs minimap2:"
python3 detect_false_positives.py \
    --nal results/adaptive_3tier.sam \
    --minimap2 "$MINIMAP2_BASELINE" \
    --tolerance 1000 \
    --nal-name "NAL-3tier"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Comparison complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Results saved to:"
echo "  • results/adaptive_2tier.sam"
echo "  • results/adaptive_3tier.sam"
echo ""

