#!/bin/bash
# Quick start script for new droid session

echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║              GenoCache V4 - New Session Quick Start                        ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""

cd /home/nebius/genocache/genocache-v4

echo "1. Verifying critical files..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

FILES=(
    "extend_phase.py"
    "adaptive_seeding.py"
    "fast_alignment.py"
    "test_extend_mock.py"
    "NEW_DROID_SESSION_HANDOFF.md"
    "HOW_NEURALIGNER_SOLVED_IT.md"
)

ALL_PRESENT=true
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✅ $file"
    else
        echo "  ❌ $file MISSING!"
        ALL_PRESENT=false
    fi
done

echo ""

if [ "$ALL_PRESENT" = false ]; then
    echo "⚠️  WARNING: Some files are missing!"
    echo "    Check CRITICAL_FILES_LIST.txt for full list"
    echo ""
fi

echo "2. Running mock test validation..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python3 test_extend_mock.py

echo ""
echo "3. Session Information"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Previous session: EXTEND phase implementation"
echo "  Status: Core fix complete ✅, validation pending ⏳"
echo "  Bug: 37% chromosome accuracy (fixed!)"
echo "  Fix: EXTEND phase (align to each, pick by score)"
echo "  Mock test: 0% → 100% accuracy ✅"
echo ""
echo "4. Next Steps"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  📖 READ: NEW_DROID_SESSION_HANDOFF.md (complete instructions)"
echo "  🔧 TODO: Test on real data (need PyTorch environment)"
echo "  🎯 GOAL: Validate 37% → 95%+ accuracy on GIAB data"
echo ""
echo "Quick commands:"
echo "  cat NEW_DROID_SESSION_HANDOFF.md    # Read full handoff"
echo "  cat HOW_NEURALIGNER_SOLVED_IT.md    # Understand the fix"
echo "  python3 test_extend_mock.py         # Re-run mock test"
echo "  python3 compare_chromosome_accuracy.py  # Check existing bug"
echo ""
echo "Ready to continue! 🚀"
echo ""

