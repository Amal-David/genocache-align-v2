#!/usr/bin/env python3
"""
Comprehensive testing of GenoCache V4 adaptive pipeline
Tests on larger sample and compares with minimap2
"""

import torch
import pickle
import faiss
from Bio import SeqIO
import sys
import time
import json
import subprocess
from pathlib import Path
from collections import defaultdict

sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

from models.encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder

print("=" * 80)
print("GenoCache V4 - Comprehensive Testing")
print("=" * 80)
print()

# Load model
print("Loading model...")
model = GenoCacheEncoder(emb_dim=128, seed_len=512)
checkpoint = torch.load('models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt', map_location='cuda')
model.load_state_dict(checkpoint['model_state_dict'])
model = model.cuda()
model.eval()
print(f"✅ Model loaded: {sum(p.numel() for p in model.parameters()):,} parameters")

# Load index  
print("Loading index...")
index = faiss.read_index('indexes/genocache_v4_production.index')
print("✅ Index loaded")

# Load metadata
print("Loading metadata...")
with open('indexes/genocache_v4_production.metadata.pkl', 'rb') as f:
    metadata = pickle.load(f)
print(f"✅ Metadata loaded: {len(metadata['positions']):,} positions")

# Initialize seeder
print("\nInitializing adaptive seeder...")
seeder = AdaptiveSeeder(model, index, metadata)
print("✅ Adaptive seeder ready")
print()

# Load test reads
print("Loading test reads...")
reads = []
for record in SeqIO.parse('/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa', 'fasta'):
    reads.append((record.id, str(record.seq)))
    if len(reads) >= 100:  # Test on 100 reads
        break
print(f"✅ Loaded {len(reads)} reads")
print()

# Run adaptive seeding
print("=" * 80)
print("ADAPTIVE SEEDING TEST - 100 Reads")
print("=" * 80)
print()

results = []
unmapped = 0
seed_counts = defaultdict(int)
status_counts = defaultdict(int)

start_time = time.time()

for read_id, read_seq in reads:
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
print("RESULTS")
print("=" * 80)
print()

print(f"Total reads: {len(reads)}")
print(f"Mapped: {len(results)} ({100*len(results)/len(reads):.1f}%)")
print(f"Unmapped: {unmapped} ({100*unmapped/len(reads):.1f}%)")
print(f"Time: {elapsed:.2f}s")
print(f"Speed: {len(reads)/elapsed:.1f} reads/sec")
print()

print("Seed usage distribution:")
for num_seeds in sorted(seed_counts.keys()):
    count = seed_counts[num_seeds]
    print(f"  {num_seeds:2d} seeds: {count:3d} reads ({100*count/len(results):.1f}%)")
print()

print("Status distribution:")
for status, count in status_counts.items():
    print(f"  {status}: {count} ({100*count/len(results):.1f}%)")
print()

# Statistics
if results:
    scores = [r['score'] for r in results]
    avg_score = sum(scores) / len(scores)
    avg_seeds = sum(r['num_seeds'] for r in results) / len(results)
    
    print("Quality metrics:")
    print(f"  Average score: {avg_score:.3f}")
    print(f"  Average seeds: {avg_seeds:.1f}")
    print(f"  Score range: {min(scores):.3f} - {max(scores):.3f}")
print()

# Save results
output = {
    'test_date': '2025-11-14',
    'num_reads': len(reads),
    'mapped': len(results),
    'unmapped': unmapped,
    'mapping_rate': len(results) / len(reads),
    'time_seconds': elapsed,
    'reads_per_second': len(reads) / elapsed,
    'seed_distribution': dict(seed_counts),
    'status_distribution': dict(status_counts),
    'avg_seeds': sum(r['num_seeds'] for r in results) / len(results) if results else 0,
    'avg_score': sum(r['score'] for r in results) / len(results) if results else 0
}

with open('adaptive_test_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print("Results saved to: adaptive_test_results.json")
print()

print("=" * 80)
print("✅ COMPREHENSIVE TEST COMPLETE!")
print("=" * 80)
print()

# Summary for minimap2 comparison
print("Key findings:")
print(f"  • Mapping rate: {100*len(results)/len(reads):.1f}%")
print(f"  • Speed: {len(reads)/elapsed:.1f} reads/sec")
print(f"  • Adaptive behavior: {min(seed_counts.keys())} - {max(seed_counts.keys())} seeds")
print(f"  • Primary alignments: {status_counts.get('primary', 0)} ({100*status_counts.get('primary', 0)/len(results):.1f}%)")
print()

# Adaptive behavior analysis
easy_reads = sum(count for seeds, count in seed_counts.items() if seeds <= 5)
hard_reads = sum(count for seeds, count in seed_counts.items() if seeds >= 16)
print(f"Adaptive strategy effectiveness:")
print(f"  • Easy reads (≤5 seeds): {easy_reads} ({100*easy_reads/len(results):.1f}%)")
print(f"  • Hard reads (≥16 seeds, rescue): {hard_reads} ({100*hard_reads/len(results):.1f}%)")
print()

print("🎯 Ready for minimap2 comparison!")
print("   Next: Run minimap2 on same reads and compare results")
