#!/bin/bash
# Start NAL index building

cd /home/nebius/genocache/genocache-v4.1-production/development/training/nal_aligned

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                      NAL INDEX BUILDING - STARTING                           ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration (following NAL paper Section A.3):"
echo "  • Seed length: 512bp (NAL uses 256-512)"
echo "  • Stride: 32bp (≤ seed_len/8 as per NAL)"
echo "  • Index type: IVFPQ (Inverted File + Product Quantization)"
echo "  • Compression: PQ16×8 (128D → 16 bytes)"
echo "  • nlist: sqrt(num_vectors) (auto-calculated)"
echo "  • nprobe: 8 (NAL recommends 8-32)"
echo "  • Distance: Inner product"
echo "  • Encoder: Use output directly (NO projection head!)"
echo ""
echo "Estimated:"
echo "  • Vectors: ~90M (GRCh38 with stride=32)"
echo "  • Memory: ~1.5GB (compressed)"
echo "  • Time: ~45 minutes on H100 GPU"
echo ""
echo "Starting in 3 seconds..."
sleep 3

source /home/nebius/genocache/.venv/bin/activate

python3 build_index_nal.py \
    --genome /home/nebius/genocache/GRCh38.fa \
    --model models/genocache_nal.pt \
    --output indexes/genocache_nal_stride32.index \
    --seed-len 512 \
    --stride 32 \
    --nprobe 8 \
    --batch-size 2048 \
    --device cuda \
    2>&1 | tee logs/indexing_$(date +%Y%m%d_%H%M%S).log

echo ""
echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                      INDEX BUILDING COMPLETE!                                ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"

