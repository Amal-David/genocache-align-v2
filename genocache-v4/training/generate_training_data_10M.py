#!/usr/bin/env python3
"""
Generate 10M training examples with improved hard negatives
"""

import sys
import h5py
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO

sys.path.append(str(Path(__file__).parent.parent))
from training.augmentation import DataAugmentation

def generate_improved_training_data(
    genome_path,
    output_path,
    num_examples=10_000_000,
    seed_len=512,
    min_negative_distance=100_000  # 100kb minimum distance for hard negatives
):
    """
    Generate training data with improved hard negatives
    
    Improvements:
    - Negatives at least 100kb away (vs 10kb before)
    - More diverse sampling across genome
    - Better quality filtering
    """
    
    print("="*80)
    print("Generating 10M Training Examples (Improved)")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Output: {output_path}")
    print(f"  Num examples: {num_examples:,}")
    print(f"  Seed length: {seed_len}")
    print(f"  Min negative distance: {min_negative_distance:,} bp")
    print()
    
    # Load genome
    print("Loading genome...")
    chromosomes = []
    chr_names = []
    chr_offsets = [0]
    
    for record in SeqIO.parse(genome_path, "fasta"):
        if record.id.startswith('NC_') or record.id.startswith('chr'):
            seq = str(record.seq).upper()
            chromosomes.append(seq)
            chr_names.append(record.id)
            chr_offsets.append(chr_offsets[-1] + len(seq))
            print(f"  Loaded {record.id}: {len(seq):,} bp")
    
    total_genome_size = chr_offsets[-1]
    print(f"\n✅ Loaded {len(chromosomes)} chromosomes")
    print(f"   Total size: {total_genome_size:,} bp\n")
    
    # Initialize augmentation
    augmenter = DataAugmentation(seed_len=seed_len)
    
    # Storage for examples
    anchors = []
    positives = []
    negatives = []
    
    print("Generating training examples...")
    valid_count = 0
    attempt_count = 0
    
    with tqdm(total=num_examples, desc="Generating") as pbar:
        while valid_count < num_examples:
            attempt_count += 1
            
            # Sample random position for anchor
            anchor_pos = np.random.randint(0, total_genome_size - seed_len)
            
            # Find which chromosome this is in
            chr_idx = 0
            for i in range(len(chr_offsets) - 1):
                if anchor_pos >= chr_offsets[i] and anchor_pos < chr_offsets[i+1]:
                    chr_idx = i
                    break
            
            local_pos = anchor_pos - chr_offsets[chr_idx]
            if local_pos + seed_len > len(chromosomes[chr_idx]):
                continue
            
            anchor_seq = chromosomes[chr_idx][local_pos:local_pos + seed_len]
            
            # Quality filter: skip if too many N's
            if anchor_seq.count('N') > seed_len * 0.1:
                continue
            
            # Generate positive with random error rate (1-10%)
            error_rate = np.random.uniform(0.01, 0.10)
            positive_seq, _ = augmenter.augment_sequence(anchor_seq, error_rate=error_rate)
            
            # Generate hard negative (at least 100kb away)
            # Try to sample from same chromosome first, then different if needed
            max_attempts = 10
            negative_found = False
            
            for _ in range(max_attempts):
                # 70% same chromosome, 30% different chromosome
                if np.random.random() < 0.7 and len(chromosomes[chr_idx]) > seed_len + min_negative_distance:
                    # Same chromosome, far away
                    neg_chr_idx = chr_idx
                    chr_len = len(chromosomes[neg_chr_idx])
                    
                    # Find position at least min_negative_distance away
                    if np.random.random() < 0.5:
                        # Try before anchor
                        max_before = max(0, local_pos - min_negative_distance - seed_len)
                        if max_before > 0:
                            neg_local_pos = np.random.randint(0, max_before)
                            negative_found = True
                    else:
                        # Try after anchor
                        min_after = local_pos + seed_len + min_negative_distance
                        if min_after + seed_len < chr_len:
                            neg_local_pos = np.random.randint(min_after, chr_len - seed_len)
                            negative_found = True
                else:
                    # Different chromosome
                    neg_chr_idx = np.random.choice([i for i in range(len(chromosomes)) if i != chr_idx])
                    chr_len = len(chromosomes[neg_chr_idx])
                    if chr_len > seed_len:
                        neg_local_pos = np.random.randint(0, chr_len - seed_len)
                        negative_found = True
                
                if negative_found:
                    break
            
            if not negative_found:
                continue
            
            negative_seq = chromosomes[neg_chr_idx][neg_local_pos:neg_local_pos + seed_len]
            
            # Quality filter for negative
            if negative_seq.count('N') > seed_len * 0.1:
                continue
            
            # Store examples
            anchors.append(anchor_seq.encode('utf-8'))
            positives.append(positive_seq.encode('utf-8'))
            negatives.append(negative_seq.encode('utf-8'))
            
            valid_count += 1
            pbar.update(1)
            
            # Periodic save (every 100k)
            if valid_count % 100_000 == 0:
                print(f"\n  Progress: {valid_count:,} / {num_examples:,} ({valid_count/num_examples*100:.1f}%)")
                print(f"  Attempt rate: {valid_count/attempt_count*100:.1f}% valid")
    
    print(f"\n✅ Generated {valid_count:,} examples")
    print(f"   Attempt rate: {valid_count/attempt_count*100:.1f}% valid\n")
    
    # Save to HDF5
    print(f"Saving to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with h5py.File(output_path, 'w') as f:
        f.create_dataset('anchors', data=anchors, dtype=h5py.string_dtype())
        f.create_dataset('positives', data=positives, dtype=h5py.string_dtype())
        f.create_dataset('negatives', data=negatives, dtype=h5py.string_dtype())
        
        # Metadata
        f.attrs['num_examples'] = num_examples
        f.attrs['seed_len'] = seed_len
        f.attrs['min_negative_distance'] = min_negative_distance
        f.attrs['genome_size'] = total_genome_size
    
    file_size = output_path.stat().st_size / (1024**3)
    print(f"✅ Saved {num_examples:,} examples ({file_size:.2f} GB)\n")
    
    print("="*80)
    print("✅ Data generation complete!")
    print("="*80)

if __name__ == "__main__":
    genome_path = Path('/home/nebius/genocache/GRCh38.fa')
    output_path = Path('data/training_10M.h5')
    
    generate_improved_training_data(
        genome_path,
        output_path,
        num_examples=10_000_000,
        min_negative_distance=100_000  # 100kb for harder negatives
    )
