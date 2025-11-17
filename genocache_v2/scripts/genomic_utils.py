"""
Genomic utilities for GenoCache-Align V2
Handles DNA sequence encoding, reference genome parsing, and data loading
"""

import torch
import numpy as np
from Bio import SeqIO
from pathlib import Path
import json
from datetime import datetime


# DNA encoding maps
DNA_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}  # N treated as A
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def one_hot_encode(sequence):
    """
    Convert DNA sequence to one-hot encoding
    
    Args:
        sequence: String of nucleotides (A, C, G, T, N)
    Returns:
        tensor: (4, len) one-hot encoded sequence
    """
    seq_len = len(sequence)
    encoding = torch.zeros(4, seq_len)
    
    for i, nucleotide in enumerate(sequence.upper()):
        idx = DNA_TO_INT.get(nucleotide, 0)
        encoding[idx, i] = 1.0
    
    return encoding


def reverse_complement(sequence):
    """Get reverse complement of DNA sequence"""
    return ''.join(COMPLEMENT[base] for base in reversed(sequence.upper()))


def extract_seeds(sequence, seed_len=512, stride=32, both_strands=False):
    """
    Extract seeds from a sequence
    
    Args:
        sequence: DNA string
        seed_len: Length of each seed
        stride: Distance between seed starts
        both_strands: If True, also extract from reverse complement
    
    Returns:
        seeds: List of (seed_seq, position, strand) tuples
    """
    seeds = []
    seq_len = len(sequence)
    
    # Forward strand
    for pos in range(0, seq_len - seed_len + 1, stride):
        seed_seq = sequence[pos:pos + seed_len]
        if len(seed_seq) == seed_len:
            seeds.append((seed_seq, pos, '+'))
    
    # Reverse strand
    if both_strands:
        rev_seq = reverse_complement(sequence)
        for pos in range(0, seq_len - seed_len + 1, stride):
            seed_seq = rev_seq[pos:pos + seed_len]
            if len(seed_seq) == seed_len:
                # Position on forward strand where this reverse seed maps
                rev_pos = seq_len - pos - seed_len
                seeds.append((seed_seq, rev_pos, '-'))
    
    return seeds


class ReferenceGenome:
    """
    Handle reference genome loading and seed extraction
    """
    
    def __init__(self, fasta_path):
        self.fasta_path = Path(fasta_path)
        self.sequences = {}
        self.load_sequences()
    
    def load_sequences(self):
        """Load all sequences from FASTA"""
        print(f"Loading reference from {self.fasta_path}...")
        
        for record in SeqIO.parse(self.fasta_path, "fasta"):
            self.sequences[record.id] = str(record.seq)
        
        total_len = sum(len(seq) for seq in self.sequences.values())
        print(f"✅ Loaded {len(self.sequences)} sequences, {total_len:,} bp total")
    
    def extract_seeds(self, chrom, seed_len=512, stride=32, both_strands=False):
        """Extract seeds from a specific chromosome"""
        if chrom not in self.sequences:
            raise ValueError(f"Chromosome {chrom} not found in reference")
        
        sequence = self.sequences[chrom]
        seeds = extract_seeds(sequence, seed_len, stride, both_strands)
        
        return seeds
    
    def get_sequence(self, chrom, start, end):
        """Extract sequence from genomic coordinates"""
        if chrom not in self.sequences:
            raise ValueError(f"Chromosome {chrom} not found in reference")
        
        return self.sequences[chrom][start:end]


class SeedDataset(torch.utils.data.Dataset):
    """
    PyTorch Dataset for seed training
    
    Returns:
        seed_tensor: (4, seed_len) one-hot encoded
        position: Genomic position (used for contrastive loss)
    """
    
    def __init__(self, seeds, positions):
        """
        Args:
            seeds: List of DNA strings
            positions: List of genomic positions
        """
        self.seeds = seeds
        self.positions = positions
    
    def __len__(self):
        return len(self.seeds)
    
    def __getitem__(self, idx):
        seed_seq = self.seeds[idx]
        position = self.positions[idx]
        
        # One-hot encode
        seed_tensor = one_hot_encode(seed_seq)
        
        return seed_tensor, position


def create_training_data(reference_path, chrom, seed_len=512, stride=32, 
                        output_dir='data/training'):
    """
    Create training dataset from reference genome
    
    Args:
        reference_path: Path to reference FASTA
        chrom: Chromosome to extract (e.g., 'NC_000022.11')
        seed_len: Seed length
        stride: Distance between seeds
        output_dir: Where to save training data
    
    Returns:
        Dataset path and metadata
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*80}")
    print("CREATING TRAINING DATA")
    print(f"{'='*80}")
    
    # Load reference
    ref = ReferenceGenome(reference_path)
    
    # Extract seeds
    print(f"\nExtracting seeds from {chrom}...")
    seeds_data = ref.extract_seeds(chrom, seed_len, stride, both_strands=False)
    
    seeds = [s[0] for s in seeds_data]
    positions = [s[1] for s in seeds_data]
    
    print(f"✅ Extracted {len(seeds):,} seeds")
    
    # Save to disk
    data_path = output_dir / f"seeds_{chrom}.npz"
    np.savez_compressed(
        data_path,
        seeds=np.array(seeds),
        positions=np.array(positions)
    )
    
    # Save metadata
    metadata = {
        'reference_path': str(reference_path),
        'chromosome': chrom,
        'seed_length': seed_len,
        'stride': stride,
        'num_seeds': len(seeds),
        'both_strands': False,
        'genomic_span': f"0-{max(positions) + seed_len}",
        'created_at': datetime.now().isoformat()
    }
    
    metadata_path = output_dir / f"seeds_{chrom}.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n✅ Saved training data:")
    print(f"   Data: {data_path}")
    print(f"   Metadata: {metadata_path}")
    
    return data_path, metadata


def load_training_data(data_path):
    """Load training data from disk"""
    data = np.load(data_path, allow_pickle=True)
    seeds = data['seeds'].tolist()
    positions = data['positions']
    
    return seeds, positions


if __name__ == "__main__":
    # Test utilities
    print("Testing genomic utilities...")
    
    # Test one-hot encoding
    seq = "ACGTACGT"
    encoded = one_hot_encode(seq)
    print(f"\nOne-hot encoding test:")
    print(f"  Input: {seq}")
    print(f"  Output shape: {encoded.shape}")
    print(f"  Sum per position: {encoded.sum(dim=0).tolist()}")  # Should be all 1s
    
    # Test reverse complement
    rev = reverse_complement(seq)
    print(f"\nReverse complement test:")
    print(f"  Input: {seq}")
    print(f"  Output: {rev}")
    
    print("\n✅ Utilities ready!")
