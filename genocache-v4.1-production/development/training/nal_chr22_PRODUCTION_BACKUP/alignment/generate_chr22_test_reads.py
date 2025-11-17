#!/usr/bin/env python3
"""
Generate synthetic Chr22 reads for testing
"""

import random

def load_chr22(fasta_path, chr_name='NC_000022.11'):
    """Load Chr22"""
    seq = []
    in_chr = False
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if chr_name in line:
                    in_chr = True
                else:
                    if in_chr:
                        break
                    in_chr = False
            elif in_chr:
                seq.append(line.upper())
    
    return ''.join(seq)

def add_ont_errors(seq, error_rate=0.10):
    """Add ONT errors"""
    result = []
    i = 0
    while i < len(seq):
        if random.random() < error_rate:
            error_type = random.choice(['sub', 'ins', 'del'])
            if error_type == 'sub':
                result.append(random.choice('ACGT'))
                i += 1
            elif error_type == 'ins':
                result.append(random.choice('ACGT'))
            else:
                i += 1
        else:
            result.append(seq[i])
            i += 1
    return ''.join(result)

def main():
    print("Generating 200 synthetic Chr22 reads...")
    
    chr_seq = load_chr22('/home/nebius/genocache/GRCh38.fa')
    print(f"Loaded Chr22: {len(chr_seq):,} bp")
    
    reads = []
    ground_truth = []
    
    for i in range(200):
        read_len = 2000
        start = random.randint(1000, len(chr_seq) - read_len - 1000)
        
        clean_seq = chr_seq[start:start + read_len]
        noisy_seq = add_ont_errors(clean_seq, 0.10)
        
        read_id = f"chr22_synthetic_{i:04d}"
        reads.append((read_id, noisy_seq))
        ground_truth.append(f"{read_id},NC_000022.11,{start},{start+read_len},{read_len},{len(noisy_seq)}\n")
        
        if (i + 1) % 50 == 0:
            print(f"  Generated {i+1}/200...")
    
    # Write FASTQ
    with open('data/chr22_test_reads.fastq', 'w') as f:
        for read_id, seq in reads:
            f.write(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")
    
    # Write ground truth
    with open('data/chr22_ground_truth.csv', 'w') as f:
        f.write("read_id,chr,start,end,clean_length,noisy_length\n")
        f.writelines(ground_truth)
    
    print("✅ Generated 200 Chr22 synthetic reads")

if __name__ == '__main__':
    main()
