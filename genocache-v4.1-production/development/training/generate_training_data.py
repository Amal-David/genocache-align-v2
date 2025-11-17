#!/usr/bin/env python3
"""
Generate training data for GenoCache V4

Strategy:
1. Load GRCh38 genome
2. Sample 130M positions (with stride)
3. Apply augmentation for each
4. Save to disk in batches (memory efficient)

Output: HDF5 file with anchors, positives, negatives
"""

import sys
import h5py
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from training.augmentation import DataAugmentation


def load_genome(genome_path, exclude_patterns=['chrM', 'chrEBV', '_random', '_alt', '_fix', 'mitochondrion']):
    """
    Load genome chromosomes (exclude mitochondria, patches, alternates)
    
    Args:
        genome_path: Path to FASTA file
        exclude_patterns: Patterns to exclude from chromosome names
    
    Returns:
        chromosomes: Dict of {name: sequence}
    """
    print(f"Loading genome from {genome_path}...")
    chromosomes = {}
    
    for record in SeqIO.parse(genome_path, "fasta"):
        # Skip unwanted chromosomes
        if any(pattern in record.id for pattern in exclude_patterns):
            continue
        
        # Keep main chromosomes only (NCBI format: NC_000001.11 or UCSC: chr1)
        if record.id.startswith('NC_') or record.id.startswith('chr'):
            seq = str(record.seq).upper()
            chromosomes[record.id] = seq
            print(f"  Loaded {record.id}: {len(seq):,} bp")
    
    total_bp = sum(len(seq) for seq in chromosomes.values())
    print(f"\n✅ Loaded {len(chromosomes)} chromosomes, {total_bp:,} bp total\n")
    
    return chromosomes


def concatenate_genome(chromosomes):
    """
    Concatenate all chromosomes into single sequence
    Add 10kb of N's between chromosomes to prevent cross-chromosome seeds
    
    Args:
        chromosomes: Dict of {name: sequence}
    
    Returns:
        genome: Single concatenated sequence
        chromosome_boundaries: List of (start, end, name) tuples
    """
    print("Concatenating genome...")
    parts = []
    boundaries = []
    current_pos = 0
    separator = 'N' * 10000
    
    for name, seq in sorted(chromosomes.items()):
        parts.append(seq)
        end_pos = current_pos + len(seq)
        boundaries.append((current_pos, end_pos, name))
        current_pos = end_pos
        
        # Add separator
        parts.append(separator)
        current_pos += len(separator)
    
    genome = ''.join(parts)
    print(f"✅ Concatenated genome: {len(genome):,} bp\n")
    
    return genome, boundaries


def generate_batch(genome, augmenter, batch_size, seed_len):
    """
    Generate one batch of training triplets
    
    Args:
        genome: Concatenated genome sequence
        augmenter: DataAugmentation instance
        batch_size: Number of triplets
        seed_len: Seed length
    
    Returns:
        anchors, positives, negatives: Training triplets
    """
    return augmenter.create_training_batch(
        genome,
        batch_size=batch_size,
        hard_negative_regions=None,  # TODO: Load V3 hard negatives
        hard_neg_ratio=0.4
    )


def save_batch_hdf5(h5file, batch_idx, anchors, positives, negatives):
    """
    Save batch to HDF5 file
    
    Args:
        h5file: h5py.File handle
        batch_idx: Batch number
        anchors, positives, negatives: Lists of sequences
    """
    # Convert to numpy arrays (as strings)
    anchors_arr = np.array(anchors, dtype='S512')  # S512 = 512-char string
    positives_arr = np.array(positives, dtype='S512')
    negatives_arr = np.array(negatives, dtype='S512')
    
    # Create datasets if first batch
    if batch_idx == 0:
        h5file.create_dataset('anchors', data=anchors_arr, maxshape=(None,), 
                             chunks=True, compression='gzip')
        h5file.create_dataset('positives', data=positives_arr, maxshape=(None,),
                             chunks=True, compression='gzip')
        h5file.create_dataset('negatives', data=negatives_arr, maxshape=(None,),
                             chunks=True, compression='gzip')
    else:
        # Append to existing datasets
        for name, arr in [('anchors', anchors_arr), 
                          ('positives', positives_arr),
                          ('negatives', negatives_arr)]:
            dataset = h5file[name]
            old_size = dataset.shape[0]
            dataset.resize(old_size + len(arr), axis=0)
            dataset[old_size:] = arr


def generate_training_data(
    genome_path,
    output_path,
    total_examples=1_000_000,  # Start with 1M for testing, scale to 130M later
    batch_size=10_000,
    seed_len=512
):
    """
    Generate full training dataset
    
    Args:
        genome_path: Path to GRCh38.fa
        output_path: Output HDF5 file
        total_examples: Total training examples to generate
        batch_size: Batch size for generation
        seed_len: Seed length
    """
    print("="*80)
    print("GenoCache V4 - Training Data Generation")
    print("="*80)
    print(f"Output: {output_path}")
    print(f"Total examples: {total_examples:,}")
    print(f"Batch size: {batch_size:,}")
    print(f"Seed length: {seed_len}")
    print("="*80)
    print()
    
    # Load genome
    chromosomes = load_genome(genome_path)
    genome, boundaries = concatenate_genome(chromosomes)
    
    # Initialize augmenter
    augmenter = DataAugmentation(seed_len=seed_len)
    
    # Generate data in batches
    num_batches = (total_examples + batch_size - 1) // batch_size
    
    print(f"Generating {num_batches} batches...\n")
    
    with h5py.File(output_path, 'w') as h5file:
        # Store metadata
        h5file.attrs['genome_path'] = str(genome_path)
        h5file.attrs['total_examples'] = total_examples
        h5file.attrs['seed_len'] = seed_len
        h5file.attrs['genome_size'] = len(genome)
        h5file.attrs['num_chromosomes'] = len(chromosomes)
        
        for batch_idx in tqdm(range(num_batches), desc="Generating batches"):
            # Generate batch
            anchors, positives, negatives = generate_batch(
                genome, augmenter, batch_size, seed_len
            )
            
            # Save to HDF5
            save_batch_hdf5(h5file, batch_idx, anchors, positives, negatives)
            
            # Progress update every 10 batches
            if (batch_idx + 1) % 10 == 0:
                total_generated = (batch_idx + 1) * batch_size
                print(f"  Generated {total_generated:,} / {total_examples:,} examples")
    
    print(f"\n✅ Saved training data to {output_path}")
    print(f"   File size: {Path(output_path).stat().st_size / 1024**2:.1f} MB")


def validate_dataset(dataset_path):
    """
    Validate generated dataset
    
    Args:
        dataset_path: Path to HDF5 file
    """
    print("\n" + "="*80)
    print("Validating dataset...")
    print("="*80)
    
    with h5py.File(dataset_path, 'r') as h5file:
        print(f"\nMetadata:")
        for key, value in h5file.attrs.items():
            print(f"  {key}: {value}")
        
        print(f"\nDatasets:")
        for name in h5file.keys():
            dataset = h5file[name]
            print(f"  {name}: shape={dataset.shape}, dtype={dataset.dtype}")
        
        # Sample a few examples
        print(f"\nSample examples:")
        anchors = h5file['anchors']
        positives = h5file['positives']
        negatives = h5file['negatives']
        
        for i in range(min(3, len(anchors))):
            print(f"\n  Example {i}:")
            print(f"    Anchor:   {anchors[i][:50]}...")
            print(f"    Positive: {positives[i][:50]}...")
            print(f"    Negative: {negatives[i][:50]}...")
    
    print("\n✅ Dataset validation complete")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate GenoCache V4 training data")
    parser.add_argument('--genome', type=str, 
                       default='/home/nebius/genocache/GRCh38.fa',
                       help='Path to GRCh38.fa')
    parser.add_argument('--output', type=str,
                       default='data/training_1M.h5',
                       help='Output HDF5 file')
    parser.add_argument('--num-examples', type=int,
                       default=1_000_000,
                       help='Number of training examples')
    parser.add_argument('--batch-size', type=int,
                       default=10_000,
                       help='Batch size')
    parser.add_argument('--seed-len', type=int,
                       default=512,
                       help='Seed length')
    parser.add_argument('--validate', action='store_true',
                       help='Validate dataset after generation')
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Generate data
    generate_training_data(
        genome_path=args.genome,
        output_path=args.output,
        total_examples=args.num_examples,
        batch_size=args.batch_size,
        seed_len=args.seed_len
    )
    
    # Validate if requested
    if args.validate:
        validate_dataset(args.output)
    
    print("\n" + "="*80)
    print("✅ Training data generation complete!")
    print("="*80)
