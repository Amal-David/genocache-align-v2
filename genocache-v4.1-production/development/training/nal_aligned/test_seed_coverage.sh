#!/bin/bash
# Test different seed numbers to find bottleneck

set -e

cd /home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned
source /home/nebius/genocache/.venv/bin/activate

BASE_DIR="/home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned"
MODEL="${BASE_DIR}/models/genocache_nal.pt"
INDEX="${BASE_DIR}/indexes/genocache_nal_stride32.index"
POSITIONS="${BASE_DIR}/indexes/genocache_nal_stride32.positions.npz"
REFERENCE="/home/nebius/genocache/GRCh38.fa"
READS="/home/nebius/genocache/genocache-v4.1-production/validation/data/giab_hg002_100reads.fastq"

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                TESTING SEED COVERAGE IMPACT                                  ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Testing different seed numbers to identify bottleneck:"
echo "  • Current: 7→13 seeds (71% mapped)"
echo "  • Test 1: 16 seeds (single iteration)"
echo "  • Test 2: 32 seeds (single iteration)"
echo "  • Test 3: 64 seeds (single iteration)"
echo ""
echo "Expected: If seeding is bottleneck, more seeds = higher mapping rate"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Test with different seed numbers
for NUM_SEEDS in 16 32 64; do
    echo ""
    echo "Testing with $NUM_SEEDS seeds..."
    OUTPUT="${BASE_DIR}/results/giab_seeds_${NUM_SEEDS}.sam"
    
    # Temporarily modify align_nal.py to use fixed seed number
    python3 << EOF
import sys
sys.path.insert(0, '${BASE_DIR}')
from seeding_nal import NALSeeding
from chaining_nal import NALChaining
from align_nal import NALAligner

# Initialize
aligner = NALAligner(
    model_path='${MODEL}',
    index_path='${INDEX}',
    positions_path='${POSITIONS}',
    reference_path='${REFERENCE}',
    device='cuda'
)

# Load reads
reads = []
with open('${READS}') as f:
    while True:
        header = f.readline().strip()
        if not header:
            break
        seq = f.readline().strip()
        plus = f.readline().strip()
        qual = f.readline().strip()
        read_id = header[1:].split()[0]
        reads.append((read_id, seq, qual))

print(f"Testing with ${NUM_SEEDS} seeds (no rescue)...")

# Align with fixed seed number (no rescue)
alignments = []
mapped = 0
for i, (read_id, read_seq, read_qual) in enumerate(reads):
    if (i + 1) % 10 == 0:
        print(f"  Progress: {i+1}/{len(reads)} reads...")
    
    read_len = len(read_seq)
    
    # Use fixed number of seeds (no rescue)
    anchors = aligner.seeder.get_anchors(read_seq, num_seeds=${NUM_SEEDS}, K=32)
    chains, _ = aligner.chainer.chain_with_rescue_check(
        anchors, read_len=read_len, num_seeds=${NUM_SEEDS}, seed_len=512
    )
    
    if chains:
        best_chain = chains[0]
        mapq = aligner._calculate_mapq(chains, best_chain)
        alignment = aligner._align_with_wfa(read_seq, best_chain, read_len)
        
        alignments.append({
            'read_id': read_id,
            'mapped': True,
            'ref_chr': best_chain['ref_chr'],
            'ref_pos': best_chain['ref_pos'],
            'strand': best_chain['strand'],
            'mapq': mapq,
            'cigar': alignment['cigar']
        })
        mapped += 1
    else:
        alignments.append({
            'read_id': read_id,
            'mapped': False
        })

print(f"Results: {mapped}/{len(reads)} = {mapped/len(reads)*100:.1f}% mapped")

# Write SAM
aligner.write_sam(alignments, '${OUTPUT}')

EOF

    # Count mapped
    TOTAL=$(grep -v "^@" "${OUTPUT}" | wc -l)
    MAPPED=$(grep -v "^@" "${OUTPUT}" | awk '$3 != "*"' | wc -l)
    RATE=$(echo "scale=1; $MAPPED * 100 / $TOTAL" | bc)
    
    echo "  → ${NUM_SEEDS} seeds: $MAPPED/$TOTAL = $RATE% mapped"
    
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "SUMMARY:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Baseline (7→13 rescue): 71%"
for NUM_SEEDS in 16 32 64; do
    OUTPUT="${BASE_DIR}/results/giab_seeds_${NUM_SEEDS}.sam"
    if [ -f "${OUTPUT}" ]; then
        TOTAL=$(grep -v "^@" "${OUTPUT}" | wc -l)
        MAPPED=$(grep -v "^@" "${OUTPUT}" | awk '$3 != "*"' | wc -l)
        RATE=$(echo "scale=1; $MAPPED * 100 / $TOTAL" | bc)
        echo "${NUM_SEEDS} seeds (no rescue): $RATE%"
    fi
done
echo ""
echo "minimap2 baseline: 94.5%"
echo ""

