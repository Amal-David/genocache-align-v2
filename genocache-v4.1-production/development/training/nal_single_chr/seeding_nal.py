#!/usr/bin/env python3
"""
NAL-Aligned Seeding - Exact NeuralAligner Protocol

Key points from NAL paper (Section 3.3):
- Seed length: 256-512bp (we use 512bp)
- Extract seeds from query read
- Encode with NAL model (128D, no projection!)
- Query FAISS index with K=32 neighbors
- Return anchors with similarity scores

Reference: NeuralAligner paper Section 3.3
"""

import sys
import torch
import faiss
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict

# Add paths
sys.path.insert(0, str(Path(__file__).parent))
from encoder_nal import NALEncoder


class NALSeeding:
    """
    NAL-aligned seeding following paper Section 3.3
    
    Uses longer seeds (512bp) with neural encoding for
    robust retrieval despite sequencing errors
    """
    
    def __init__(
        self,
        model_path: str,
        index_path: str,
        positions_path: str,
        seed_len: int = 512,
        K: int = 32,  # NAL uses K=32 neighbors
        device: str = 'cuda'
    ):
        """
        Initialize NAL seeding
        
        Args:
            model_path: Path to trained NAL model
            index_path: Path to FAISS index
            positions_path: Path to positions file
            seed_len: Seed length (NAL: 256-512)
            K: Number of nearest neighbors (NAL: 32)
            device: Device for model
        """
        self.seed_len = seed_len
        self.K = K
        self.device = device
        
        # Load model
        print(f"Loading NAL model from {model_path}...")
        checkpoint = torch.load(model_path, map_location=device)
        self.model = NALEncoder(
            emb_dim=128,
            seed_len=seed_len,
            vocab_size=5,
            hidden_dim=128,
            num_layers=4
        ).to(device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.eval()
        if 'epoch' in checkpoint:
            print(f"  Model loaded (epoch {checkpoint['epoch']}, val_loss {checkpoint['val_metrics']['loss']:.4f})")
        else:
            print(f"  Model loaded (batch {checkpoint.get('batch', 'unknown')}, loss {checkpoint.get('loss', 0):.4f})")
        
        # Load index
        print(f"Loading FAISS index from {index_path}...")
        self.index = faiss.read_index(index_path)
        print(f"  Index loaded: {self.index.ntotal:,} vectors")
        
        # Load positions
        print(f"Loading positions from {positions_path}...")
        pos_data = np.load(positions_path, allow_pickle=True)
        self.chr_names = pos_data['chr_names']
        self.positions = pos_data['positions']
        print(f"  Positions loaded: {len(self.positions):,} entries")
        
        # Base to index mapping
        self.base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def extract_seeds_from_read(
        self,
        read_seq: str,
        num_seeds: int = 5,
        strategy: str = 'uniform'
    ) -> List[Tuple[int, str]]:
        """
        Extract seeds from query read (NAL Section 3.3, 3.5)
        
        NAL adaptive strategy:
        - First iteration: 5 seeds
        - Rescue iteration: 16 seeds (if needed)
        
        Args:
            read_seq: Query read sequence
            num_seeds: Number of seeds to extract
            strategy: 'uniform' or 'random'
        
        Returns:
            List of (position, seed_seq) tuples
        """
        read_len = len(read_seq)
        
        if read_len < self.seed_len:
            # Read too short, pad
            seed = read_seq + 'N' * (self.seed_len - read_len)
            return [(0, seed)]
        
        seeds = []
        
        if strategy == 'uniform':
            # Uniformly spaced seeds
            if num_seeds == 1:
                pos = (read_len - self.seed_len) // 2
                seeds.append((pos, read_seq[pos:pos + self.seed_len]))
            else:
                step = max(1, (read_len - self.seed_len) // (num_seeds - 1))
                for i in range(num_seeds):
                    pos = min(i * step, read_len - self.seed_len)
                    seeds.append((pos, read_seq[pos:pos + self.seed_len]))
        
        else:  # random
            import random
            positions = sorted(random.sample(
                range(read_len - self.seed_len + 1),
                min(num_seeds, read_len - self.seed_len + 1)
            ))
            for pos in positions:
                seeds.append((pos, read_seq[pos:pos + self.seed_len]))
        
        return seeds
    
    def encode_seeds(self, seed_seqs: List[str]) -> np.ndarray:
        """
        Encode seeds to embeddings (NAL: use encoder output, no projection!)
        
        CRITICAL: use_projection=False (Section A.3)
        
        Args:
            seed_seqs: List of DNA sequences
        
        Returns:
            Embeddings (N, 128)
        """
        # Convert to tensors
        batch_tensors = []
        for seq in seed_seqs:
            indices = [self.base_to_idx.get(b.upper(), 4) for b in seq[:self.seed_len]]
            # Pad if needed
            if len(indices) < self.seed_len:
                indices = indices + [4] * (self.seed_len - len(indices))
            batch_tensors.append(torch.tensor(indices, dtype=torch.long))
        
        x = torch.stack(batch_tensors).to(self.device)
        
        # Encode (no projection head!)
        with torch.no_grad():
            embeddings = self.model(x, use_projection=False)
        
        return embeddings.cpu().numpy()
    
    def query_index(
        self,
        embeddings: np.ndarray,
        K: int = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Query FAISS index for nearest neighbors (NAL Section 3.3)
        
        Args:
            embeddings: Query embeddings (N, 128)
            K: Number of neighbors (default: self.K)
        
        Returns:
            distances: (N, K) similarity scores
            indices: (N, K) vector indices
        """
        if K is None:
            K = self.K
        
        # Normalize for inner product (NAL uses inner product)
        embeddings = embeddings.astype('float32')
        faiss.normalize_L2(embeddings)
        
        # Search
        distances, indices = self.index.search(embeddings, K)
        
        return distances, indices
    
    def get_anchors(
        self,
        read_seq: str,
        num_seeds: int = 5,
        K: int = None
    ) -> List[Dict]:
        """
        Complete seeding: extract → encode → search → anchors
        
        Args:
            read_seq: Query read sequence
            num_seeds: Number of seeds to extract
            K: Number of neighbors per seed
        
        Returns:
            List of anchor dicts with:
                - seed_idx: Seed index in read
                - seed_pos: Position of seed in read
                - ref_chr: Reference chromosome
                - ref_pos: Reference position
                - similarity: Similarity score
                - anchor_idx: Index in results (for tracking)
        """
        if K is None:
            K = self.K
        
        # Extract seeds
        seeds = self.extract_seeds_from_read(read_seq, num_seeds)
        
        if not seeds:
            return []
        
        # Encode seeds
        seed_seqs = [s[1] for s in seeds]
        embeddings = self.encode_seeds(seed_seqs)
        
        # Query index
        distances, indices = self.query_index(embeddings, K)
        
        # Convert to anchors
        anchors = []
        for seed_idx, (seed_pos, seed_seq) in enumerate(seeds):
            for anchor_idx in range(K):
                vec_idx = indices[seed_idx, anchor_idx]
                similarity = distances[seed_idx, anchor_idx]
                
                # Skip invalid indices
                if vec_idx < 0 or vec_idx >= len(self.positions):
                    continue
                
                # Get genomic position
                chr_name = self.chr_names[vec_idx]
                ref_pos = self.positions[vec_idx]
                
                anchors.append({
                    'seed_idx': seed_idx,
                    'seed_pos': seed_pos,
                    'ref_chr': chr_name,
                    'ref_pos': int(ref_pos),
                    'similarity': float(similarity),
                    'anchor_idx': anchor_idx
                })
        
        return anchors


def test_seeding():
    """Test NAL seeding"""
    import sys
    
    if len(sys.argv) < 4:
        print("Usage: python seeding_nal.py <model> <index> <positions>")
        print("\nTest example:")
        print("  python seeding_nal.py \\")
        print("    models/genocache_nal.pt \\")
        print("    indexes/genocache_nal_stride32.index \\")
        print("    indexes/genocache_nal_stride32.positions.npz")
        return
    
    model_path = sys.argv[1]
    index_path = sys.argv[2]
    positions_path = sys.argv[3]
    
    print("="* 80)
    print("NAL Seeding Test")
    print("=" * 80)
    
    # Initialize
    seeder = NALSeeding(
        model_path=model_path,
        index_path=index_path,
        positions_path=positions_path,
        seed_len=512,
        K=32
    )
    
    # Test read (mock)
    test_read = "ACGT" * 300  # 1200bp read
    
    print(f"\nTest read: {len(test_read)}bp")
    print(f"Extracting 5 seeds...")
    
    # Get anchors
    anchors = seeder.get_anchors(test_read, num_seeds=5, K=32)
    
    print(f"\nResults:")
    print(f"  Seeds extracted: 5")
    print(f"  Anchors found: {len(anchors)}")
    print(f"  Anchors per seed: {len(anchors) / 5:.1f}")
    
    # Show top anchor for each seed
    print(f"\nTop anchor per seed:")
    for seed_idx in range(5):
        seed_anchors = [a for a in anchors if a['seed_idx'] == seed_idx]
        if seed_anchors:
            top = max(seed_anchors, key=lambda x: x['similarity'])
            print(f"  Seed {seed_idx}: {top['ref_chr']}:{top['ref_pos']} "
                  f"(sim={top['similarity']:.3f})")
    
    print("\n✅ Seeding test complete!")


if __name__ == '__main__':
    test_seeding()
