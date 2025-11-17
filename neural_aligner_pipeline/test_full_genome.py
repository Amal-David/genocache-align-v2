#!/usr/bin/env python3
"""
Phase 5.4: Test full genome alignment pipeline
Validates the complete system with synthetic reads
"""

import sys
import os
import numpy as np
from Bio import SeqIO
import random
import time

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from multi_seeder import MultiSeedAligner
from refine_alignment import RefinedAligner
from quality_scorer import QualityScorer
from bam_writer import BAMWriter

# Configuration
REF_FA = "/home/nebius/genocache/GRCh38.fa"
ENCODER_PT = "/home/nebius/genocache/nal_encoder_full_genome.pt"
FAISS_INDEX = "/home/nebius/genocache/grch38_full_index.faiss"
POSITIONS_NPY = "/home/nebius/genocache/grch38_full_positions.npy"
OUTPUT_BAM = "/home/nebius/genocache/test_full_genome_aligned.bam"

# Test parameters
NUM_TEST_READS = 100  # Start with 100 reads for quick test
READ_LENGTH = 1000
ERROR_RATE = 0.05

print("="*70)
print("Phase 5.4: Full Genome Alignment Test")
print("="*70)

# Step 1: Load reference
print("\n[1/6] Loading reference genome...")
ref_seqs = {}
chr_lengths = {}

for record in SeqIO.parse(REF_FA, "fasta"):
    chr_name = record.id
    if not chr_name.startswith('NC_'):
        continue
    
    seq = str(record.seq).upper()
    ref_seqs[chr_name] = seq
    chr_lengths[chr_name] = len(seq)
    print(f"  ✓ Loaded {chr_name}: {len(seq):,} bp")

print(f"\n  Total chromosomes: {len(ref_seqs)}")
print(f"  Total genome size: {sum(chr_lengths.values()) / 1e9:.2f} Gbp")

# Step 2: Generate synthetic test reads
print(f"\n[2/6] Generating {NUM_TEST_READS} synthetic test reads...")

test_reads = []
ground_truth = []

for i in range(NUM_TEST_READS):
    # Pick random chromosome (weighted by length)
    total_length = sum(chr_lengths.values())
    chr_probs = {chr: length/total_length for chr, length in chr_lengths.items()}
    chr_name = random.choices(list(chr_probs.keys()), weights=list(chr_probs.values()))[0]
    
    # Pick random position
    chr_seq = ref_seqs[chr_name]
    if len(chr_seq) < READ_LENGTH + 1000:
        continue
    
    true_pos = random.randint(1000, len(chr_seq) - READ_LENGTH - 1000)
    
    # Extract read
    read_seq = chr_seq[true_pos:true_pos + READ_LENGTH]
    
    # Add sequencing errors
    read_seq = list(read_seq)
    bases = ['A', 'C', 'G', 'T']
    for j in range(len(read_seq)):
        if random.random() < ERROR_RATE:
            # Substitution
            read_seq[j] = random.choice([b for b in bases if b != read_seq[j]])
    
    read_seq = ''.join(read_seq)
    
    # Skip if too many Ns
    if read_seq.count('N') > READ_LENGTH * 0.1:
        continue
    
    read_name = f"read_{i}_{chr_name}_{true_pos}"
    test_reads.append((read_name, read_seq))
    ground_truth.append((chr_name, true_pos))

print(f"  ✓ Generated {len(test_reads)} valid test reads")

# Step 3: Initialize alignment pipeline
print("\n[3/6] Initializing full genome alignment pipeline...")

# Load multi-seeder
print("  Loading multi-seed aligner...")
multi_seeder = MultiSeedAligner(
    encoder_path=ENCODER_PT,
    index_path=FAISS_INDEX,
    positions_path=POSITIONS_NPY
)
print("  ✓ Multi-seeder ready")

# Load SW refiner (simple)
print("  Loading SW refiner...")
from refine_alignment import AlignmentRefiner
sw_refiner = AlignmentRefiner(ref_fa=REF_FA)
print("  ✓ SW refiner ready")

# Quality scorer
scorer = QualityScorer()
print("  ✓ Quality scorer ready")

# BAM writer (pass reference file path, not dict)
bam_writer = BAMWriter(OUTPUT_BAM, REF_FA)
print("  ✓ BAM writer ready")

# Step 4: Align all reads
print(f"\n[4/6] Aligning {len(test_reads)} reads...")

results = []
start_time = time.time()

for i, ((read_name, read_seq), (true_chr, true_pos)) in enumerate(zip(test_reads, ground_truth)):
    # Multi-seed alignment
    candidates = multi_seeder.align_read(read_seq, top_k=5)
    
    if not candidates:
        results.append({
            'read_name': read_name,
            'aligned': False,
            'true_chr': true_chr,
            'true_pos': true_pos
        })
        continue
    
    # Take top candidate
    chr_name, pos, score = candidates[0]
    
    # Refine with Smith-Waterman
    refined = sw_refiner.refine_alignment(read_seq, pos, chromosome=chr_name)
    
    if refined:
        refined_pos = refined['position']
        identity = refined['identity']
        cigar = refined['cigar']
        
        # Calculate MAPQ
        mapq = scorer.calculate_mapq(identity, score, len(candidates))
        
        # Write to BAM
        bam_writer.write_alignment(
            query_name=read_name,
            query_sequence=read_seq,
            reference_name=chr_name,
            reference_start=refined_pos,
            cigar_string=cigar,
            mapping_quality=mapq
        )
        
        # Calculate error
        error = abs(refined_pos - true_pos) if chr_name == true_chr else 999999
        
        results.append({
            'read_name': read_name,
            'aligned': True,
            'true_chr': true_chr,
            'true_pos': true_pos,
            'pred_chr': chr_name,
            'pred_pos': refined_pos,
            'error_bp': error,
            'mapq': mapq,
            'identity': identity
        })
    else:
        results.append({
            'read_name': read_name,
            'aligned': False,
            'true_chr': true_chr,
            'true_pos': true_pos
        })
    
    if (i + 1) % 10 == 0:
        elapsed = time.time() - start_time
        rate = (i + 1) / elapsed
        print(f"  Progress: {i+1}/{len(test_reads)} ({rate:.1f} reads/sec)")

elapsed = time.time() - start_time
bam_writer.close()

print(f"  ✓ Aligned in {elapsed:.1f}s ({len(test_reads)/elapsed:.1f} reads/sec)")

# Step 5: Calculate metrics
print("\n[5/6] Calculating accuracy metrics...")

aligned = [r for r in results if r['aligned']]
correct_chr = [r for r in aligned if r['pred_chr'] == r['true_chr']]
errors = [r['error_bp'] for r in correct_chr]

recall = len(aligned) / len(results) * 100
chr_accuracy = len(correct_chr) / len(aligned) * 100 if aligned else 0

if errors:
    median_error = np.median(errors)
    exact_matches = sum(1 for e in errors if e == 0)
    within_100bp = sum(1 for e in errors if e <= 100)
    within_1kb = sum(1 for e in errors if e <= 1000)
    
    exact_pct = exact_matches / len(errors) * 100
    within_100_pct = within_100bp / len(errors) * 100
    within_1kb_pct = within_1kb / len(errors) * 100
else:
    median_error = 0
    exact_pct = 0
    within_100_pct = 0
    within_1kb_pct = 0

# Step 6: Print results
print("\n[6/6] Results Summary")
print("="*70)
print(f"\nTest Configuration:")
print(f"  Test reads:        {len(test_reads)}")
print(f"  Read length:       {READ_LENGTH} bp")
print(f"  Error rate:        {ERROR_RATE*100:.1f}%")
print(f"  Genome:            Full GRCh38 ({len(ref_seqs)} chromosomes)")

print(f"\nAlignment Results:")
print(f"  Aligned:           {len(aligned)}/{len(results)} ({recall:.1f}%)")
print(f"  Correct chr:       {len(correct_chr)}/{len(aligned)} ({chr_accuracy:.1f}%)")

print(f"\nPosition Accuracy (correct chromosome only):")
print(f"  Median error:      {median_error:.0f} bp")
print(f"  Exact matches:     {exact_pct:.1f}%")
print(f"  Within 100bp:      {within_100_pct:.1f}%")
print(f"  Within 1kb:        {within_1kb_pct:.1f}%")

if aligned:
    avg_mapq = np.mean([r['mapq'] for r in aligned])
    avg_identity = np.mean([r['identity'] for r in aligned])
    print(f"\nQuality Metrics:")
    print(f"  Average MAPQ:      {avg_mapq:.1f}")
    print(f"  Average identity:  {avg_identity:.1f}%")

print(f"\nPerformance:")
print(f"  Total time:        {elapsed:.1f}s")
print(f"  Throughput:        {len(test_reads)/elapsed:.1f} reads/sec")

print(f"\nOutput:")
print(f"  BAM file:          {OUTPUT_BAM}")

print("\n" + "="*70)
print("✓ Phase 5.4 Test Complete!")
print("="*70)

# Save results
results_file = "/home/nebius/genocache/test_full_genome_results.txt"
with open(results_file, 'w') as f:
    f.write("Full Genome Alignment Test Results\n")
    f.write("="*70 + "\n\n")
    f.write(f"Recall: {recall:.1f}%\n")
    f.write(f"Chromosome accuracy: {chr_accuracy:.1f}%\n")
    f.write(f"Median error: {median_error:.0f} bp\n")
    f.write(f"Exact matches: {exact_pct:.1f}%\n")
    f.write(f"Within 100bp: {within_100_pct:.1f}%\n")
    f.write(f"Throughput: {len(test_reads)/elapsed:.1f} reads/sec\n")

print(f"\nResults saved to: {results_file}")
