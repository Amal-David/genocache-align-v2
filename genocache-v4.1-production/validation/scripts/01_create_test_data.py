#!/usr/bin/env python3
"""
Create test dataset from chr22 reference

Generates 100 realistic reads from chr22 for validation
"""

import random
from Bio import SeqIO
from pathlib import Path

def create_test_reads(reference_path, output_path, num_reads=100, read_length=1000):
    """
    Create test reads by extracting regions from reference
    Add realistic sequencing errors (~2%)
    """
    print("="*80)
    print("Creating Test Dataset from chr22 Reference")
    print("="*80)
    print()
    
    # Load chr22
    print(f"Loading reference from {reference_path}...")
    sequences = {}
    for record in SeqIO.parse(reference_path, "fasta"):
        if 'chr22' in record.id.lower() or record.id.startswith('NC_000022'):
            sequences[record.id] = str(record.seq).upper()
            print(f"  Loaded {record.id}: {len(record.seq):,} bp")
    
    if not sequences:
        print("❌ No chr22 found!")
        return
    
    # Get main chr22
    chr22_id = list(sequences.keys())[0]
    chr22_seq = sequences[chr22_id]
    
    print(f"\n✅ Using {chr22_id} ({len(chr22_seq):,} bp)")
    print(f"\nGenerating {num_reads} reads...")
    
    # Generate reads
    reads = []
    ground_truth = []
    random.seed(42)  # Reproducible
    
    for i in range(num_reads):
        # Skip N-rich regions
        attempts = 0
        while attempts < 100:
            # Random position
            start = random.randint(10000, len(chr22_seq) - read_length - 10000)
            end = start + read_length
            
            read_seq = chr22_seq[start:end]
            
            # Check quality
            if read_seq.count('N') < read_length * 0.05:
                break
            attempts += 1
        
        if attempts >= 100:
            continue
        
        # Add sequencing errors (~2% error rate)
        read_list = list(read_seq)
        num_errors = int(read_length * 0.02)
        
        for _ in range(num_errors):
            pos = random.randint(0, len(read_list) - 1)
            error_type = random.choice(['sub', 'sub', 'sub', 'ins', 'del'])  # More subs
            
            if error_type == 'sub':
                read_list[pos] = random.choice(['A', 'C', 'G', 'T'])
            elif error_type == 'ins' and len(read_list) < read_length * 1.1:
                read_list.insert(pos, random.choice(['A', 'C', 'G', 'T']))
            elif error_type == 'del' and len(read_list) > read_length * 0.9:
                if pos < len(read_list):
                    del read_list[pos]
        
        read_seq_with_errors = ''.join(read_list)
        
        # Store
        read_id = f"test_read_{i:03d}"
        reads.append((read_id, read_seq_with_errors))
        ground_truth.append({
            'read_id': read_id,
            'chr': chr22_id,
            'start': start,
            'end': end,
            'length': len(read_seq_with_errors)
        })
        
        if (i + 1) % 20 == 0:
            print(f"  Generated {i + 1}/{num_reads} reads...")
    
    print(f"\n✅ Generated {len(reads)} reads")
    
    # Write FASTQ
    output_path = Path(output_path)
    print(f"\nWriting FASTQ to {output_path}...")
    with open(output_path, 'w') as f:
        for read_id, seq in reads:
            # FASTQ format
            f.write(f"@{read_id}\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{'I' * len(seq)}\n")  # Quality scores (all high quality)
    
    print(f"✅ Wrote {len(reads)} reads")
    
    # Write ground truth
    truth_path = str(output_path).replace('.fastq', '_ground_truth.txt')
    print(f"\nWriting ground truth to {truth_path}...")
    with open(truth_path, 'w') as f:
        f.write("read_id\tchr\tstart\tend\tlength\n")
        for gt in ground_truth:
            f.write(f"{gt['read_id']}\t{gt['chr']}\t{gt['start']}\t{gt['end']}\t{gt['length']}\n")
    
    print(f"✅ Wrote ground truth")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total reads: {len(reads)}")
    print(f"Read length: ~{read_length} bp (with errors)")
    print(f"Error rate: ~2%")
    print(f"Chromosome: {chr22_id}")
    print()
    print(f"Output files:")
    print(f"  Reads: {output_path}")
    print(f"  Truth: {truth_path}")
    
    return reads, ground_truth


if __name__ == '__main__':
    # Paths
    reference_path = Path('/home/nebius/genocache/GRCh38.fa')
    output_path = Path(__file__).parent.parent / 'data' / 'test_reads_100.fastq'
    
    if not reference_path.exists():
        print(f"❌ Reference not found: {reference_path}")
        exit(1)
    
    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate test data
    reads, truth = create_test_reads(reference_path, output_path, num_reads=100, read_length=1000)
    
    print("\n✅ Test dataset ready for validation!")
