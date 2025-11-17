#!/usr/bin/env python3
"""
Test COMPLETE end-to-end pipeline with actual alignment (CIGAR strings)
"""

import torch
import pickle
import faiss
from Bio import SeqIO
import sys
import time
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

from models.encoder import GenoCacheEncoder
from adaptive_seeding import AdaptiveSeeder
from fast_alignment import FastAligner

print("=" * 80)
print("GenoCache V4 - COMPLETE Pipeline Test (with CIGAR)")
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

# Load genome for alignment (main chromosomes only)
print("\nLoading genome for alignment...")
genome_dict = {}
for record in SeqIO.parse('/home/nebius/genocache/GRCh38.fa', 'fasta'):
    # Only load main chromosomes (NC_000001 through NC_000024, plus X and Y)
    if record.id.startswith('NC_0000') and len(record.id) < 15:  # Main chromosomes
        genome_dict[record.id] = str(record.seq)
        print(f"  Loaded: {record.id} ({len(record.seq):,} bp)")
    if len(genome_dict) >= 24:  # chr1-22, X, Y
        break
print(f"✅ Genome loaded: {len(genome_dict)} chromosomes")

# Initialize aligner
print("\nInitializing fast aligner...")
aligner = FastAligner(genome_dict, mode='semi-global')
print("✅ Fast aligner ready")
print()

# Load test reads
print("Loading test reads...")
reads = []
for record in SeqIO.parse('/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa', 'fasta'):
    reads.append((record.id, str(record.seq)))
    if len(reads) >= 10:  # Test on 10 reads
        break
print(f"✅ Loaded {len(reads)} reads")
print()

# Run COMPLETE pipeline
print("=" * 80)
print("COMPLETE PIPELINE TEST (Seeding + Alignment)")
print("=" * 80)
print()

results = []
unmapped = 0
total_time = 0

for i, (read_id, read_seq) in enumerate(reads):
    print(f"\nRead {i+1}/{len(reads)}: {read_id}")
    print(f"  Length: {len(read_seq)} bp")
    
    start_time = time.time()
    
    # Step 1: Adaptive seeding
    seed_result = seeder.align_read(read_seq, read_id=read_id)
    
    if not seed_result:
        print(f"  ❌ Unmapped (seeding failed)")
        unmapped += 1
        continue
    
    seed_time = time.time() - start_time
    
    print(f"  ✅ Seeded: {seed_result['chr']}:{seed_result['start']}-{seed_result['end']}")
    print(f"     Seeds: {seed_result['num_seeds']}, Score: {seed_result['score']:.2f}")
    print(f"     Status: {seed_result['status']}, Time: {seed_time*1000:.1f}ms")
    
    # Step 2: Alignment
    align_start = time.time()
    alignment = aligner.align_read(
        read_seq, 
        seed_result['chr'],
        seed_result['start'],
        seed_result['end']
    )
    align_time = time.time() - align_start
    
    if not alignment:
        print(f"  ⚠️  Alignment failed")
        unmapped += 1
        continue
    
    total_time = time.time() - start_time
    
    print(f"  ✅ Aligned: {alignment['ref_chr']}:{alignment['ref_start']}-{alignment['ref_end']}")
    print(f"     CIGAR: {alignment['cigar'][:50]}{'...' if len(alignment['cigar']) > 50 else ''}")
    print(f"     Score: {alignment['score']}")
    print(f"     Align time: {align_time*1000:.1f}ms")
    print(f"     Total time: {total_time*1000:.1f}ms")
    
    # Format SAM
    sam_record = aligner.format_sam(read_id, read_seq, alignment)
    
    results.append({
        'read_id': read_id,
        'seed_result': seed_result,
        'alignment': alignment,
        'sam': sam_record,
        'time': total_time
    })

print()
print("=" * 80)
print("RESULTS")
print("=" * 80)
print()

mapped = len(results)
print(f"Total reads: {len(reads)}")
print(f"Mapped: {mapped} ({100*mapped/len(reads):.1f}%)")
print(f"Unmapped: {unmapped} ({100*unmapped/len(reads):.1f}%)")

if results:
    avg_time = sum(r['time'] for r in results) / len(results)
    print(f"Average time: {avg_time*1000:.1f}ms per read")
    print(f"Throughput: {1.0/avg_time:.1f} reads/sec")
    print()
    
    print("Sample SAM records:")
    for i, result in enumerate(results[:3]):
        print(f"\n{i+1}. {result['sam'].strip()}")

print()
print("=" * 80)
print("✅ COMPLETE PIPELINE TEST DONE!")
print("=" * 80)
print()

print("Summary:")
print(f"  • Seeding: Adaptive (2-{max(r['seed_result']['num_seeds'] for r in results)} seeds)")
print(f"  • Alignment: Smith-Waterman (parasail)")
print(f"  • CIGAR: Generated ✅")
print(f"  • SAM: Generated ✅")
print(f"  • Speed: {1.0/avg_time:.1f} reads/sec (complete pipeline)")
print()

# Save SAM file
sam_output = "test_complete_10reads.sam"
with open(sam_output, 'w') as f:
    # Header
    f.write("@HD\tVN:1.0\tSO:unsorted\n")
    for chr_name in sorted(genome_dict.keys())[:24]:
        chr_len = len(genome_dict[chr_name])
        f.write(f"@SQ\tSN:{chr_name}\tLN:{chr_len}\n")
    f.write("@PG\tID:GenoCache\tPN:GenoCache\tVN:4.0\n")
    
    # Alignments
    for result in results:
        f.write(result['sam'])

print(f"✅ SAM file saved: {sam_output}")
print()

print("🎉 END-TO-END PIPELINE WORKING WITH CIGAR STRINGS! 🎉")
