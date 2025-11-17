"""
GenoCache V4 - Dataset with Augmentation
Implements error simulation, position shifts, and reverse complement
"""

import torch
from torch.utils.data import Dataset
import numpy as np
from Bio import SeqIO
import random

class DNAEncoder:
    """Encode DNA sequences to integers"""
    def __init__(self):
        self.char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        self.idx_to_char = {v: k for k, v in self.char_to_idx.items()}
    
    def encode(self, seq):
        """Convert DNA string to list of integers"""
        return [self.char_to_idx.get(c.upper(), 4) for c in seq]
    
    def decode(self, indices):
        """Convert integers back to DNA string"""
        return ''.join(self.idx_to_char.get(i, 'N') for i in indices)

class GenomeDataset(Dataset):
    """
    Genome dataset with augmentation for contrastive learning
    
    Augmentations:
    - Error simulation (1-10%): substitutions, insertions, deletions
    - Position shift (±51bp): to learn translation invariance
    - Reverse complement (50%): to learn strand invariance
    """
    
    def __init__(self, 
                 genome_path, 
                 window_size=512, 
                 stride=256,
                 max_windows=None,
                 error_rate_range=(0.01, 0.10),
                 position_shift=51,
                 rc_prob=0.5,
                 chromosomes=None,
                 min_n_content=0.9):
        """
        Args:
            genome_path: Path to genome FASTA file
            window_size: Size of each window in bp (512)
            stride: Stride between windows (256 for training, 32 for indexing)
            max_windows: Maximum number of windows to generate (None = all)
            error_rate_range: Range of error rates (0.01, 0.10) = 1-10%
            position_shift: Max position shift in bp (±51)
            rc_prob: Probability of reverse complement (0.5)
            chromosomes: List of chromosome names to include (None = all primary)
            min_n_content: Minimum non-N content (0.9 = skip if >10% N)
        """
        self.window_size = window_size
        self.stride = stride
        self.error_rate_range = error_rate_range
        self.position_shift = position_shift
        self.rc_prob = rc_prob
        self.min_n_content = min_n_content
        self.encoder = DNAEncoder()
        
        print(f"\nLoading genome from {genome_path}...")
        self.sequences = {}
        self.chr_order = []
        
        for record in SeqIO.parse(genome_path, "fasta"):
            # Filter chromosomes if specified
            if chromosomes is not None:
                if record.id not in chromosomes:
                    continue
            
            # Skip non-primary chromosomes
            if any(x in record.id.lower() for x in 
                   ['alt', 'random', 'un', 'fix', 'patch', 'hap']):
                continue
            
            self.sequences[record.id] = str(record.seq).upper()
            self.chr_order.append(record.id)
            print(f"  Loaded {record.id}: {len(record.seq):,} bp")
        
        print(f"\nTotal chromosomes: {len(self.sequences)}")
        print(f"Total genome size: {sum(len(s) for s in self.sequences.values()):,} bp")
        
        # Generate window positions
        print(f"\nGenerating windows (window={window_size}, stride={stride})...")
        self.windows = []
        
        for chr_name in self.chr_order:
            seq = self.sequences[chr_name]
            chr_windows = 0
            
            for pos in range(0, len(seq) - window_size, stride):
                # Check N content
                window = seq[pos:pos + window_size]
                non_n_ratio = sum(1 for c in window if c != 'N') / len(window)
                
                if non_n_ratio < self.min_n_content:
                    continue
                
                self.windows.append((chr_name, pos))
                chr_windows += 1
                
                # Check max windows limit
                if max_windows and len(self.windows) >= max_windows:
                    break
            
            print(f"  {chr_name}: {chr_windows:,} windows")
            
            if max_windows and len(self.windows) >= max_windows:
                break
        
        print(f"\nTotal windows: {len(self.windows):,}")
        print(f"Memory reduction: {(len(seq) / len(self.windows)):.1f}× (stride={stride})")
    
    def __len__(self):
        return len(self.windows)
    
    def augment_sequence(self, seq):
        """
        Apply error augmentation
        
        Simulates ONT sequencing errors:
        - Substitutions (33%)
        - Insertions (33%)
        - Deletions (33%)
        """
        seq_list = list(seq)
        error_rate = random.uniform(*self.error_rate_range)
        num_errors = int(len(seq) * error_rate)
        
        for _ in range(num_errors):
            if len(seq_list) == 0:
                break
            
            pos = random.randint(0, len(seq_list) - 1)
            error_type = random.random()
            
            if error_type < 0.33:  # Substitution
                seq_list[pos] = random.choice(['A', 'C', 'G', 'T'])
            elif error_type < 0.66:  # Deletion
                if len(seq_list) > 1:
                    seq_list.pop(pos)
            else:  # Insertion
                seq_list.insert(pos, random.choice(['A', 'C', 'G', 'T']))
        
        # Ensure exactly window_size
        result = ''.join(seq_list)
        if len(result) > self.window_size:
            result = result[:self.window_size]
        elif len(result) < self.window_size:
            result = result + 'N' * (self.window_size - len(result))
        
        return result
    
    def reverse_complement(self, seq):
        """Compute reverse complement"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(c, 'N') for c in reversed(seq))
    
    def __getitem__(self, idx):
        """
        Get one training example
        
        Returns:
            anchor: Augmented DNA sequence (encoded as integers)
            chr_name: Chromosome name (for debugging)
            pos: Position in genome (for debugging)
        """
        chr_name, pos = self.windows[idx]
        seq = self.sequences[chr_name]
        
        # Apply random position shift
        shift = random.randint(-self.position_shift, self.position_shift)
        start = max(0, pos + shift)
        end = min(len(seq), start + self.window_size)
        
        # Extract window
        window = seq[start:end]
        
        # Pad if needed (at chromosome boundaries)
        if len(window) < self.window_size:
            window = window + 'N' * (self.window_size - len(window))
        
        # Apply error augmentation
        window_aug = self.augment_sequence(window)
        
        # Apply reverse complement with probability
        if random.random() < self.rc_prob:
            window_aug = self.reverse_complement(window_aug)
        
        # Encode to integers
        anchor_enc = torch.tensor(self.encoder.encode(window_aug), dtype=torch.long)
        
        return anchor_enc, chr_name, pos

def collate_fn(batch):
    """Custom collate function for DataLoader"""
    anchors, chr_names, positions = zip(*batch)
    anchors = torch.stack(anchors)
    return anchors, chr_names, positions

if __name__ == '__main__':
    # Test dataset
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python dataset.py <genome.fa>")
        sys.exit(1)
    
    genome_path = sys.argv[1]
    
    print("Creating dataset...")
    dataset = GenomeDataset(
        genome_path=genome_path,
        window_size=512,
        stride=256,
        max_windows=1000,
        error_rate_range=(0.01, 0.10),
        position_shift=51,
        rc_prob=0.5
    )
    
    print(f"\nDataset size: {len(dataset)}")
    
    print("\nTesting getitem...")
    anchor, chr_name, pos = dataset[0]
    print(f"  Anchor shape: {anchor.shape}")
    print(f"  Chromosome: {chr_name}")
    print(f"  Position: {pos}")
    print(f"  First 20 tokens: {anchor[:20].tolist()}")
    
    print("\nTesting DataLoader...")
    from torch.utils.data import DataLoader
    loader = DataLoader(dataset, batch_size=4, shuffle=True, collate_fn=collate_fn)
    
    for anchors, chr_names, positions in loader:
        print(f"  Batch shape: {anchors.shape}")
        print(f"  Chromosomes: {chr_names}")
        print(f"  Positions: {positions}")
        break
    
    print("\n✅ Dataset test passed!")
