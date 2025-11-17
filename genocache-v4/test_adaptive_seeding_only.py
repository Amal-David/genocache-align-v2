#!/usr/bin/env python3
"""
Quick test of adaptive seeding only (skip WFA for speed)
"""

import torch
import pickle
import faiss
from Bio import SeqIO
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

from models.encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder

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
    if len(reads) >= 10:
        break
print(f"✅ Loaded {len(reads)} reads")
print()

# Test adaptive seeding
print("Testing adaptive seeding...")
print("=" * 80)

for read_id, read_seq in reads:
    print(f"\nRead: {read_id} ({len(read_seq)} bp)")
    
    result = seeder.align_read(read_seq, read_id=read_id)
    
    if result:
        print(f"  ✅ Mapped to {result['chr']}:{result['start']}-{result['end']}")
        print(f"     Seeds used: {result['num_seeds']}")
        print(f"     Score: {result['score']:.3f}")
        print(f"     Status: {result['status']}")
    else:
        print(f"  ❌ Unmapped")

print()
print("=" * 80)
print("✅ Adaptive seeding test complete!")
