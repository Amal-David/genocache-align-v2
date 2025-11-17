#!/bin/bash
# Quick start script for GenoCache v6.1 full genome training

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║              GenoCache v6.1 - Full Genome Training                           ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  • Genome: GRCh38 (all 24 chromosomes)"
echo "  • Batches: 24,000"
echo "  • Curriculum: 5% → 15% errors"
echo "  • Expected time: 4-6 hours"
echo ""
echo "Based on Chr22 success: 76.5% accuracy (beat minimap2!)"
echo ""
read -p "Start training now? (y/n) " -n 1 -r
echo ""

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Starting training..."
    cd training/
    python3 train_full_genome.py 2>&1 | tee ../logs/training_$(date +%Y%m%d_%H%M%S).log
else
    echo "Cancelled."
fi
