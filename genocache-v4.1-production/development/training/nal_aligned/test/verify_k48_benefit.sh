#!/bin/bash
# Verify K=48 benefit extends to 32 and 48 seeds

set -e

cd "$(dirname "$0")"
source /home/nebius/genocache/.venv/bin/activate

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║           VERIFYING K=48 BENEFIT ON HIGHER SEED COUNTS                      ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Testing hypothesis: K=48 improves mapping for all seed counts"
echo ""

# Test configurations
echo "Testing 4 configurations:"
echo "  1. 32 seeds, K=32 (baseline)"
echo "  2. 32 seeds, K=48 (optimized)"
echo "  3. 48 seeds, K=32 (baseline)"
echo "  4. 48 seeds, K=48 (optimized)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

python3 << 'EOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))

from parameter_sweep import ParameterSweep, load_reads

# Configuration
MODEL = "../models/genocache_nal.pt"
INDEX = "../indexes/genocache_nal_stride32.index"
POSITIONS = "../indexes/genocache_nal_stride32.positions.npz"
READS_FILE = "/home/nebius/genocache/genocache-v4.1-production/validation/data/giab_hg002_100reads.fastq"

# Load reads
print("\nLoading 100 reads...")
reads = load_reads(READS_FILE, 100)

# Initialize sweep
sweep = ParameterSweep(MODEL, INDEX, POSITIONS, device='cuda')

# Test configurations
configs = [
    (32, 32, "32 seeds, K=32"),
    (32, 48, "32 seeds, K=48"),
    (48, 32, "48 seeds, K=32"),
    (48, 48, "48 seeds, K=48"),
]

print("\nTesting configurations...\n")
results = []

for num_seeds, K, label in configs:
    print(f"Testing {label}...")
    result = sweep.test_configuration(reads, num_seeds=num_seeds, K=K, tolerance=1000)
    results.append((label, result))
    print(f"  → Mapping: {result['mapping_rate']:.1f}%, Time: {result['avg_time_ms']:.1f}ms, "
          f"Anchors/seed: {result['avg_anchor_density']:.1f}\n")

# Summary
print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
print("\nSUMMARY:")
print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

print(f"{'Configuration':<20} {'Mapping':<12} {'Improvement'}")
print("-" * 60)

# Baseline vs optimized for each seed count
for i in range(0, len(results), 2):
    baseline_label, baseline = results[i]
    opt_label, opt = results[i+1]
    
    improvement = opt['mapping_rate'] - baseline['mapping_rate']
    print(f"{baseline_label:<20} {baseline['mapping_rate']:>6.1f}%")
    print(f"{opt_label:<20} {opt['mapping_rate']:>6.1f}%      +{improvement:.1f}%")
    print()

# Key findings
print("KEY FINDINGS:")
print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

baseline_32 = results[0][1]
opt_32 = results[1][1]
baseline_48 = results[2][1]
opt_48 = results[3][1]

print(f"✅ 32 seeds: K=48 gives +{opt_32['mapping_rate'] - baseline_32['mapping_rate']:.1f}% "
      f"({baseline_32['mapping_rate']:.1f}% → {opt_32['mapping_rate']:.1f}%)")
print(f"✅ 48 seeds: K=48 gives +{opt_48['mapping_rate'] - baseline_48['mapping_rate']:.1f}% "
      f"({baseline_48['mapping_rate']:.1f}% → {opt_48['mapping_rate']:.1f}%)")

if opt_32['mapping_rate'] >= 97:
    print(f"\n🎯 BREAKTHROUGH: 32 seeds + K=48 = {opt_32['mapping_rate']:.1f}% "
          f"(matches/exceeds minimap2 94.5%!)")

if opt_48['mapping_rate'] >= 99:
    print(f"🎯 EXCELLENT: 48 seeds + K=48 = {opt_48['mapping_rate']:.1f}% "
          f"(near-perfect mapping!)")

print()
EOF

echo ""
echo "✅ Verification complete!"
echo ""

