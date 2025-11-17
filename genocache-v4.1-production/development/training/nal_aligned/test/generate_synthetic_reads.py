#!/usr/bin/env python3
"""
Generate synthetic reads from reference genome to test NAL accuracy

This creates reads with known ground truth positions,
allowing us to definitively test if NAL can map them correctly.
"""

import random
import sys
from pathlib import Path

def load_reference_chr(fasta_path, chr_name, max_len=None):
    """Load a chromosome from FASTA"""
    print(f"Loading {chr_name} from {fasta_path}...")
    seq = []
    in_chr = False
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if chr_name in line:
                    in_chr = True
                    print(f"  Found {chr_name}")
                else:
                    if in_chr:
                        break
                    in_chr = False
            elif in_chr:
                seq.append(line.upper())
                if max_len and len(''.join(seq)) >= max_len:
                    break
    
    full_seq = ''.join(seq)
    print(f"  Loaded {len(full_seq):,} bp")
    return full_seq

def add_ont_errors(seq, error_rate=0.10):
    """Add realistic ONT sequencing errors (substitutions, insertions, deletions)"""
    result = []
    i = 0
    
    while i < len(seq):
        if random.random() < error_rate:
            error_type = random.choice(['sub', 'ins', 'del'])
            
            if error_type == 'sub':
                # Substitution
                result.append(random.choice('ACGT'))
                i += 1
            elif error_type == 'ins':
                # Insertion
                result.append(random.choice('ACGT'))
                # Don't advance i (insert without consuming)
            else:  # del
                # Deletion (skip base)
                i += 1
        else:
            result.append(seq[i])
            i += 1
    
    return ''.join(result)

def generate_reads(reference_path, num_reads=200, read_length=2000):
    """Generate synthetic reads from reference"""
    print(f"\nGenerating {num_reads} synthetic reads (length ~{read_length}bp)...")
    
    # Load Chr1 (largest chromosome, lots of space)
    chr_seq = load_reference_chr(reference_path, 'NC_000001.11', max_len=50_000_000)
    
    if len(chr_seq) < read_length * 2:
        print(f"❌ Chromosome too short ({len(chr_seq)} bp)")
        return None
    
    reads = []
    ground_truth = []
    
    for i in range(num_reads):
        # Random position (avoid edges)
        max_start = len(chr_seq) - read_length - 1000
        start_pos = random.randint(1000, max_start)
        end_pos = start_pos + read_length
        
        # Extract clean sequence
        clean_seq = chr_seq[start_pos:end_pos]
        
        # Add ONT errors
        noisy_seq = add_ont_errors(clean_seq, error_rate=0.10)
        
        # Create read
        read_id = f"synthetic_read_{i:04d}"
        reads.append((read_id, noisy_seq))
        
        # Ground truth
        ground_truth.append({
            'read_id': read_id,
            'chr': 'NC_000001.11',
            'start': start_pos,
            'end': end_pos,
            'clean_length': len(clean_seq),
            'noisy_length': len(noisy_seq)
        })
        
        if (i + 1) % 50 == 0:
            print(f"  Generated {i+1}/{num_reads} reads...")
    
    return reads, ground_truth

def write_fastq(reads, output_path):
    """Write reads to FASTQ format"""
    with open(output_path, 'w') as f:
        for read_id, seq in reads:
            f.write(f"@{read_id}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * len(seq) + "\n")  # Fake quality scores

def write_ground_truth(ground_truth, output_path):
    """Write ground truth positions"""
    with open(output_path, 'w') as f:
        f.write("read_id,chr,start,end,clean_length,noisy_length\n")
        for gt in ground_truth:
            f.write(f"{gt['read_id']},{gt['chr']},{gt['start']},{gt['end']},"
                   f"{gt['clean_length']},{gt['noisy_length']}\n")

def main():
    reference_path = '/home/nebius/genocache/GRCh38.fa'
    output_fastq = 'data/synthetic_200reads.fastq'
    output_truth = 'data/synthetic_ground_truth.csv'
    
    print("="*80)
    print("GENERATING SYNTHETIC READS FOR NAL VALIDATION")
    print("="*80)
    
    # Generate reads
    result = generate_reads(reference_path, num_reads=200, read_length=2000)
    
    if result is None:
        print("❌ Failed to generate reads")
        return 1
    
    reads, ground_truth = result
    
    # Write outputs
    print(f"\nWriting FASTQ to {output_fastq}...")
    write_fastq(reads, output_fastq)
    
    print(f"Writing ground truth to {output_truth}...")
    write_ground_truth(ground_truth, output_truth)
    
    print("\n" + "="*80)
    print("✅ SYNTHETIC DATA GENERATED")
    print("="*80)
    print(f"\nCreated:")
    print(f"  • {len(reads)} synthetic reads")
    print(f"  • All from Chr1 (NC_000001.11)")
    print(f"  • ~2000bp each with ~10% ONT errors")
    print(f"  • Ground truth positions recorded")
    print(f"\nFiles:")
    print(f"  • {output_fastq}")
    print(f"  • {output_truth}")
    print(f"\nNext steps:")
    print(f"  1. Run NAL 2-tier on synthetic_200reads.fastq")
    print(f"  2. Run minimap2 on synthetic_200reads.fastq")
    print(f"  3. Compare both against ground_truth.csv")
    print(f"  4. This will show if NAL can map reads with known positions")
    print()
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
