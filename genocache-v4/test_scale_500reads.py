#!/usr/bin/env python3
"""
Scale test: Run on all 500 reads to validate at production scale
"""

import torch
import pickle
import faiss
from Bio import SeqIO
import sys
import time
import json
from pathlib import Path
from collections import defaultdict

sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

from models.encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder

print("=" * 80)
print("GenoCache V4 - SCALE TEST (500 Reads)")
print("=" * 80)
print()

# Load model
print("Loading model...")
model = GenoCacheEncoder(emb_dim=128, seed_len=512)
checkpoint = torch.load('models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt', map_location='cuda')
model.load_state_dict(checkpoint['model_state_dict'])
model = model.cuda()
model.eval()
print(f"✅ Model loaded")

# Load index  
print("Loading index...")
index = faiss.read_index('indexes/genocache_v4_production.index')
print("✅ Index loaded")

# Load metadata
print("Loading metadata...")
with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
    metadata = pickle.load(f)
print(f"✅ Metadata loaded")

# Initialize seeder
print("\nInitializing adaptive seeder...")
seeder = AdaptiveSeeder(model, index, metadata)
print("✅ Adaptive seeder ready")
print()

# Load ALL test reads
print("Loading all test reads...")
reads = []
for record in SeqIO.parse('/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa', 'fasta'):
    reads.append((record.id, str(record.seq)))
print(f"✅ Loaded {len(reads)} reads")
print()

# Run adaptive seeding
print("=" * 80)
print(f"PROCESSING {len(reads)} READS...")
print("=" * 80)
print()

results = []
unmapped = 0
seed_counts = defaultdict(int)
status_counts = defaultdict(int)

start_time = time.time()

# Process with progress updates
for i, (read_id, read_seq) in enumerate(reads):
    if (i + 1) % 50 == 0:
        elapsed = time.time() - start_time
        speed = (i + 1) / elapsed
        print(f"  Processed {i+1}/{len(reads)} reads ({speed:.1f} reads/sec)")
    
    result = seeder.align_read(read_seq, read_id=read_id)
    
    if result:
        results.append(result)
        seed_counts[result['num_seeds']] += 1
        status_counts[result['status']] += 1
    else:
        unmapped += 1

elapsed = time.time() - start_time

print()
print("=" * 80)
print("SCALE TEST RESULTS")
print("=" * 80)
print()

print(f"Total reads: {len(reads)}")
print(f"Mapped: {len(results)} ({100*len(results)/len(reads):.1f}%)")
print(f"Unmapped: {unmapped} ({100*unmapped/len(reads):.1f}%)")
print(f"Time: {elapsed:.2f}s ({elapsed/60:.1f} min)")
print(f"Speed: {len(reads)/elapsed:.1f} reads/sec")
print()

print("Seed usage statistics:")
all_seeds = [r['num_seeds'] for r in results]
print(f"  Min: {min(all_seeds)} seeds")
print(f"  Max: {max(all_seeds)} seeds")
print(f"  Average: {sum(all_seeds)/len(all_seeds):.1f} seeds")
print(f"  Median: {sorted(all_seeds)[len(all_seeds)//2]} seeds")
print()

print("Status distribution:")
for status, count in status_counts.items():
    print(f"  {status}: {count} ({100*count/len(results):.1f}%)")
print()

# Adaptive behavior analysis
easy_reads = sum(count for seeds, count in seed_counts.items() if seeds <= 5)
medium_reads = sum(count for seeds, count in seed_counts.items() if 6 <= seeds <= 15)
hard_reads = sum(count for seeds, count in seed_counts.items() if seeds >= 16)

print("Adaptive strategy breakdown:")
print(f"  Easy reads (≤5 seeds):     {easy_reads:3d} ({100*easy_reads/len(results):.1f}%)")
print(f"  Medium reads (6-15 seeds): {medium_reads:3d} ({100*medium_reads/len(results):.1f}%)")
print(f"  Hard reads (≥16 seeds):    {hard_reads:3d} ({100*hard_reads/len(results):.1f}%)")
print()

# Quality metrics
scores = [r['score'] for r in results]
print("Quality metrics:")
print(f"  Average score: {sum(scores)/len(scores):.3f}")
print(f"  Score range: {min(scores):.3f} - {max(scores):.3f}")
print()

# Save results
output = {
    'test_date': '2025-11-14',
    'test_type': 'scale_test_500_reads',
    'num_reads': len(reads),
    'mapped': len(results),
    'unmapped': unmapped,
    'mapping_rate': len(results) / len(reads),
    'time_seconds': elapsed,
    'reads_per_second': len(reads) / elapsed,
    'seed_stats': {
        'min': min(all_seeds),
        'max': max(all_seeds),
        'average': sum(all_seeds) / len(all_seeds),
        'median': sorted(all_seeds)[len(all_seeds)//2]
    },
    'status_distribution': dict(status_counts),
    'adaptive_breakdown': {
        'easy': easy_reads,
        'medium': medium_reads,
        'hard': hard_reads
    }
}

with open('scale_test_500reads_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print("=" * 80)
print("✅ SCALE TEST COMPLETE!")
print("=" * 80)
print()

print("Key findings at scale:")
print(f"  • Mapped {100*len(results)/len(reads):.1f}% of {len(reads)} reads")
print(f"  • Throughput: {len(reads)/elapsed:.1f} reads/sec")
print(f"  • Adaptive: {min(all_seeds)}-{max(all_seeds)} seeds (avg {sum(all_seeds)/len(all_seeds):.1f})")
print(f"  • Primary alignments: {status_counts.get('primary', 0)} ({100*status_counts.get('primary', 0)/len(results):.1f}%)")
print()

print("Results saved to: scale_test_500reads_results.json")
print()
print("🎯 Production-scale validation complete!")
