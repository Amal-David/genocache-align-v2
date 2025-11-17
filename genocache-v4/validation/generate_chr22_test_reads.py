#!/usr/bin/env python3
"""
Generate synthetic test reads from chr22 for validation

Creates 512bp reads with 95% identity (typical ONT)
Records ground truth positions for accuracy measurement
"""

import random
import sys
from pathlib import Path
from typing import List, Tuple
import json

sys.path.append(str(Path(__file__).parent.parent / "training"))
from augmentation import DataAugmentation


def load_chr22(genome_path: Path) -> str:
    """Load chr22 from genome FASTA"""
    print(f"Loading chr22 from {genome_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(genome_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr is not None:
                    genome[current_chr] = ''.join(current_seq).upper()
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_chr is not None:
            genome[current_chr] = ''.join(current_seq).upper()
    
    # Try different chr22 names
    for name in ['NC_000022.11', 'chr22', '22']:
        if name in genome:
            print(f"✅ Found chr22 as '{name}': {len(genome[name]):,} bp")
            return genome[name]
    
    raise ValueError("Could not find chr22 in genome!")


def generate_test_reads(
    chr22_seq: str,
    num_reads: int = 1000,
    read_length: int = 512,
    error_rate: float = 0.05,
    seed: int = 42
) -> List[Tuple[str, int, int]]:
    """
    Generate synthetic test reads
    
    Args:
        chr22_seq: chr22 sequence
        num_reads: Number of reads to generate
        read_length: Length of each read (512bp to match training)
        error_rate: Error rate (0.05 = 95% identity)
        seed: Random seed for reproducibility
    
    Returns:
        List of (read_seq, start_pos, end_pos) tuples
    """
    random.seed(seed)
    augmenter = DataAugmentation(seed_len=read_length)
    
    reads = []
    chr_len = len(chr22_seq)
    
    print(f"Generating {num_reads} test reads...")
    print(f"  Read length: {read_length}bp")
    print(f"  Error rate: {error_rate} ({100*(1-error_rate):.1f}% identity)")
    
    attempts = 0
    max_attempts = num_reads * 3
    
    while len(reads) < num_reads and attempts < max_attempts:
        attempts += 1
        
        # Sample random position
        start_pos = random.randint(0, chr_len - read_length)
        end_pos = start_pos + read_length
        
        # Extract clean sequence
        clean_seq = chr22_seq[start_pos:end_pos]
        
        # Skip if too many N's
        if clean_seq.count('N') > read_length * 0.1:
            continue
        
        # Add errors (simulating sequencing)
        noisy_seq = augmenter.add_errors(
            clean_seq,
            error_rate=error_rate,
            sub_prob=0.6,  # ONT profile
            ins_prob=0.2,
            del_prob=0.2
        )
        
        # Trim/pad to exact length
        if len(noisy_seq) > read_length:
            noisy_seq = noisy_seq[:read_length]
        elif len(noisy_seq) < read_length:
            noisy_seq = noisy_seq + 'N' * (read_length - len(noisy_seq))
        
        # 50% chance of reverse complement
        if random.random() < 0.5:
            noisy_seq = augmenter.reverse_complement(noisy_seq)
        
        reads.append((noisy_seq, start_pos, end_pos))
        
        if len(reads) % 100 == 0:
            print(f"  Generated {len(reads)}/{num_reads} reads...")
    
    print(f"✅ Generated {len(reads)} reads")
    return reads


def save_test_reads(reads: List[Tuple[str, int, int]], output_dir: Path):
    """Save test reads to files"""
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Save sequences (FASTA-like)
    seq_file = output_dir / "test_reads.fasta"
    with open(seq_file, 'w') as f:
        for i, (seq, start, end) in enumerate(reads):
            f.write(f">read_{i} pos={start}-{end}\n")
            f.write(f"{seq}\n")
    
    print(f"✅ Saved sequences: {seq_file}")
    
    # Save ground truth positions (JSON)
    truth_file = output_dir / "ground_truth.json"
    truth_data = []
    for i, (seq, start, end) in enumerate(reads):
        truth_data.append({
            'read_id': f'read_{i}',
            'true_start': start,
            'true_end': end,
            'length': len(seq)
        })
    
    with open(truth_file, 'w') as f:
        json.dump(truth_data, f, indent=2)
    
    print(f"✅ Saved ground truth: {truth_file}")
    
    # Statistics
    positions = [start for _, start, _ in reads]
    print(f"\nStatistics:")
    print(f"  Min position: {min(positions):,}")
    print(f"  Max position: {max(positions):,}")
    print(f"  Median position: {sorted(positions)[len(positions)//2]:,}")
    print(f"  Total chr22 coverage: {(max(positions) - min(positions)):,} bp")


def main():
    print("=" * 80)
    print("GenoCache V4 - Generate chr22 Test Reads")
    print("=" * 80)
    print()
    
    # Configuration
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    OUTPUT_DIR = Path("/home/nebius/genocache/genocache-v4/validation/data")
    
    NUM_READS = 1000
    READ_LENGTH = 512  # Match training!
    ERROR_RATE = 0.05  # 95% identity (typical ONT)
    
    # Load chr22
    chr22_seq = load_chr22(GENOME_PATH)
    
    # Generate reads
    reads = generate_test_reads(
        chr22_seq,
        num_reads=NUM_READS,
        read_length=READ_LENGTH,
        error_rate=ERROR_RATE
    )
    
    # Save
    save_test_reads(reads, OUTPUT_DIR)
    
    print("\n✅ Test read generation complete!")
    print(f"\nNext step: Run validation script with trained model")


if __name__ == "__main__":
    main()
