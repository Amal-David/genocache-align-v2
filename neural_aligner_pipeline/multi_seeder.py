#!/usr/bin/env python3
"""
Phase 1: Multi-seed extraction and consensus alignment for chr22
Improves 70.6% → 80-85% by using multiple seeds per read
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import faiss
import random

# Model architecture (embedded to avoid import issues)
class ImprovedNAL(nn.Module):
    def __init__(self, out_dim=256):
        super().__init__()
        self.in_proj = nn.Conv1d(4, 128, kernel_size=9, padding=4)
        self.bn1 = nn.BatchNorm1d(128)
        self.conv1 = nn.Conv1d(128, 256, kernel_size=7, padding=3)
        self.bn2 = nn.BatchNorm1d(256)
        self.conv2 = nn.Conv1d(256, 256, kernel_size=7, padding=3)
        self.bn3 = nn.BatchNorm1d(256)
        self.conv3 = nn.Conv1d(256, 512, kernel_size=5, padding=2)
        self.bn4 = nn.BatchNorm1d(512)
        self.conv4 = nn.Conv1d(512, 512, kernel_size=5, padding=2)
        self.bn5 = nn.BatchNorm1d(512)
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.proj = nn.Sequential(nn.Linear(512, out_dim), nn.LayerNorm(out_dim))
    
    def forward(self, x):
        x = x.permute(0,2,1)
        x = F.relu(self.bn1(self.in_proj(x)))
        identity = x
        x = F.relu(self.bn2(self.conv1(x)))
        x = self.bn3(self.conv2(x))
        if identity.shape[1] != x.shape[1]:
            identity = F.conv1d(identity, torch.eye(256, 128, device=x.device).unsqueeze(2), padding=0)
        x = F.relu(x + identity)
        x = F.relu(self.bn4(self.conv3(x)))
        identity = x
        x = self.bn5(self.conv4(x))
        x = F.relu(x + identity)
        x = self.pool(x).view(x.shape[0], -1)
        x = self.proj(x)
        return F.normalize(x, dim=-1)

class MultiSeedAligner:
    def __init__(self, encoder_path, index_path, positions_path, device="cuda"):
        """Initialize multi-seed aligner"""
        print(f"Loading multi-seed aligner...")
        
        # Load model
        self.device = device if torch.cuda.is_available() else "cpu"
        self.encoder = ImprovedNAL().to(self.device)
        self.encoder.load_state_dict(torch.load(encoder_path, map_location=self.device))
        self.encoder.eval()
        print(f"  ✓ Model loaded on {self.device}")
        
        # Load index
        self.index = faiss.read_index(index_path)
        print(f"  ✓ FAISS index loaded: {self.index.ntotal:,} vectors")
        
        # Load positions
        self.ref_positions = np.load(positions_path)
        print(f"  ✓ Position map loaded: {len(self.ref_positions):,} entries")
        
        self.seed_len = 256
        print(f"  ✓ Multi-seeder ready")
    
    def align_read(self, read_seq, n_seeds=5, top_k=50, min_score=0.5, cluster_dist=500):
        """
        Align read using multiple seeds
        
        Args:
            read_seq: Read sequence string
            n_seeds: Number of seeds to extract from read
            top_k: Number of top candidates per seed
            min_score: Minimum similarity score to consider
            cluster_dist: Distance threshold for clustering hits
        
        Returns:
            List of candidate alignments sorted by score
        """
        # Extract seeds from different positions
        seeds = self.extract_seeds(read_seq, n_seeds)
        
        if not seeds:
            return []
        
        # Search each seed
        all_hits = []
        
        for read_pos, seed_seq in seeds:
            # Encode seed
            vec = self.encode_seed(seed_seq)
            
            # Search index
            scores, indices = self.index.search(vec.reshape(1, -1), k=top_k)
            
            # Store hits
            for idx, score in zip(indices[0], scores[0]):
                if score < min_score:
                    continue
                
                ref_pos = self.ref_positions[idx]
                
                # Estimate read start position
                # If seed is at position X in read, and matches ref position Y,
                # then read likely starts at Y - X
                estimated_start = ref_pos - read_pos
                
                all_hits.append({
                    'ref_pos': estimated_start,
                    'seed_read_pos': read_pos,
                    'seed_ref_pos': ref_pos,
                    'score': score,
                    'raw_idx': idx
                })
        
        # Cluster nearby hits
        clusters = self.cluster_hits(all_hits, cluster_dist)
        
        return clusters
    
    def extract_seeds(self, read_seq, n_seeds, strategy='evenly_spaced'):
        """
        Extract multiple seeds from a read
        
        Strategies:
        - 'evenly_spaced': Seeds at regular intervals
        - 'high_quality': Seeds from regions with fewer Ns
        """
        read_len = len(read_seq)
        
        # Read too short
        if read_len < self.seed_len:
            return [(0, read_seq)]
        
        seeds = []
        
        if strategy == 'evenly_spaced':
            # Calculate step size
            if read_len <= self.seed_len * n_seeds:
                # Need overlapping seeds
                step = max(1, (read_len - self.seed_len) // (n_seeds - 1)) if n_seeds > 1 else 0
            else:
                # Can have non-overlapping seeds
                step = (read_len - self.seed_len) // n_seeds
            
            for i in range(n_seeds):
                pos = min(i * step, read_len - self.seed_len)
                seed = read_seq[pos:pos + self.seed_len]
                
                # Skip seeds with too many Ns
                if seed.count('N') < self.seed_len * 0.3:
                    seeds.append((pos, seed))
        
        elif strategy == 'high_quality':
            # Score all possible seeds by quality
            candidate_seeds = []
            step = self.seed_len // 2  # Overlap by 50%
            
            for i in range(0, read_len - self.seed_len + 1, step):
                seed = read_seq[i:i + self.seed_len]
                n_count = seed.count('N')
                quality = (self.seed_len - n_count) / self.seed_len
                candidate_seeds.append((i, seed, quality))
            
            # Take top N by quality
            candidate_seeds.sort(key=lambda x: x[2], reverse=True)
            seeds = [(pos, seq) for pos, seq, qual in candidate_seeds[:n_seeds]]
        
        return seeds
    
    def cluster_hits(self, hits, cluster_dist):
        """
        Group hits that are nearby (likely same alignment)
        
        Args:
            hits: List of hit dictionaries
            cluster_dist: Maximum distance to consider hits in same cluster
        
        Returns:
            List of cluster summaries sorted by total score
        """
        if not hits:
            return []
        
        # Sort by estimated read start position
        hits.sort(key=lambda h: h['ref_pos'])
        
        clusters = []
        current_cluster = [hits[0]]
        
        for hit in hits[1:]:
            last = current_cluster[-1]
            
            # Same region?
            if abs(hit['ref_pos'] - last['ref_pos']) < cluster_dist:
                current_cluster.append(hit)
            else:
                # Save current cluster and start new one
                clusters.append(self.summarize_cluster(current_cluster))
                current_cluster = [hit]
        
        # Don't forget last cluster
        if current_cluster:
            clusters.append(self.summarize_cluster(current_cluster))
        
        # Sort by total score (more seeds matching = higher confidence)
        clusters.sort(key=lambda c: c['total_score'], reverse=True)
        
        return clusters
    
    def summarize_cluster(self, hits):
        """Summarize a cluster of hits"""
        return {
            'ref_pos': int(np.median([h['ref_pos'] for h in hits])),
            'n_seeds': len(hits),
            'total_score': sum(h['score'] for h in hits),
            'avg_score': np.mean([h['score'] for h in hits]),
            'min_score': min(h['score'] for h in hits),
            'max_score': max(h['score'] for h in hits),
            'hits': hits
        }
    
    def encode_seed(self, seed_seq):
        """Encode seed sequence to vector"""
        onehot = self.seq_to_onehot(seed_seq)
        
        with torch.no_grad():
            tensor = torch.from_numpy(onehot).unsqueeze(0).to(self.device)
            vec = self.encoder(tensor).cpu().numpy()
        
        return vec
    
    def seq_to_onehot(self, seq):
        """Convert sequence to one-hot encoding"""
        arr = np.zeros((len(seq), 4), dtype=np.float32)
        
        for i, ch in enumerate(seq.upper()):
            if ch == "A":
                arr[i, 0] = 1
            elif ch == "C":
                arr[i, 1] = 1
            elif ch == "G":
                arr[i, 2] = 1
            elif ch == "T":
                arr[i, 3] = 1
            elif ch == "N":
                # Ambiguous base - use small values for all
                arr[i, random.randrange(4)] = 0.25
            else:
                # Unknown base - random
                arr[i, random.randrange(4)] = 1.0
        
        return arr


def main():
    """Quick test of multi-seeder"""
    print("="*70)
    print("Multi-Seed Aligner - Phase 1")
    print("="*70)
    
    # Create aligner
    aligner = MultiSeedAligner(
        encoder_path="nal_encoder_best.pt",
        index_path="faiss_index_improved.idx",
        positions_path="ref_positions_improved.npy"
    )
    
    # Test with synthetic read
    print(f"\nTesting with 4kb synthetic read...")
    test_read = "ACGT" * 1000
    
    results = aligner.align_read(test_read, n_seeds=5, top_k=20)
    
    print(f"\nMulti-seed alignment results:")
    print(f"  Found {len(results)} candidate alignments")
    
    for i, cluster in enumerate(results[:5]):
        print(f"  {i+1}. Position: {cluster['ref_pos']:,}")
        print(f"     Seeds matched: {cluster['n_seeds']}")
        print(f"     Total score: {cluster['total_score']:.3f}")
        print(f"     Avg score: {cluster['avg_score']:.3f}")
    
    print("\n" + "="*70)
    print("✓ Multi-seeder test complete")
    print("="*70)


if __name__ == "__main__":
    main()
