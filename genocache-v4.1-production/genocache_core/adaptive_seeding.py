"""
GenoCache V4 - Adaptive Seeding Strategy
Implements NeuralAligner-style adaptive seeding with rescue seeds
"""

import torch
import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

@dataclass
class Seed:
    """Represents a single seed from a read"""
    read_pos: int           # Position in read
    seq: str               # 512bp sequence
    embedding: np.ndarray  # 128D embedding
    candidates: List[Tuple[str, int, float]]  # (chr, pos, score)
    status: str = "pending"  # pending, unique, ambiguous, repeat

@dataclass
class SeedChain:
    """Represents a chain of seeds mapping to same chromosome"""
    chr_name: str
    seeds: List[Seed]
    score: float
    colinear: bool
    start_pos: int
    end_pos: int

class AdaptiveSeeder:
    """
    Adaptive seeding strategy from NeuralAligner
    
    Strategy:
    1. Extract 5 evenly-spaced seeds from read
    2. Search each seed in FAISS index
    3. Filter seeds (unique/ambiguous/repeat)
    4. Chain seeds using colinearity check
    5. Add rescue seeds if needed (up to 16 total)
    """
    
    def __init__(self, 
                 model,
                 index,
                 metadata,
                 min_seeds: int = 5,
                 max_seeds: int = 16,
                 window_size: int = 512,
                 top_k: int = 32,
                 uniqueness_threshold: float = 0.1,
                 colinearity_tolerance: int = 3000):
        """
        Args:
            model: Neural encoder model
            index: FAISS index
            metadata: Index metadata (positions, chr_names)
            min_seeds: Minimum number of seeds (5)
            max_seeds: Maximum number of seeds (16)
            window_size: Seed size in bp (512)
            top_k: Number of candidates per seed (32)
            uniqueness_threshold: Score difference for uniqueness
            colinearity_tolerance: Max distance error for chaining (3kb, relaxed for real reads)
        """
        self.model = model
        self.index = index
        self.metadata = metadata
        self.min_seeds = min_seeds
        self.max_seeds = max_seeds
        self.window_size = window_size
        self.top_k = top_k
        self.uniqueness_threshold = uniqueness_threshold
        self.colinearity_tolerance = colinearity_tolerance
        
        # CRITICAL FIX: Set nprobe to NeuralAligner's range (8-32)
        # Testing shows nprobe=64 is too high for real data
        # Lower values = more focused search, less noise
        if hasattr(self.index, 'nprobe'):
            original_nprobe = self.index.nprobe
            # Set to 24 for balance (NAL uses 8-32)
            self.index.nprobe = 24
            print(f"  FAISS nprobe: {original_nprobe} → 24 (optimized for real data)")
        
        # DNA encoding
        self.char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        
        self.model.eval()
    
    def encode_sequence(self, seq: str) -> np.ndarray:
        """Encode DNA sequence to embedding (128D or 256D depending on model)"""
        # Convert to indices
        indices = [self.char_to_idx.get(c.upper(), 4) for c in seq]
        
        # Pad/truncate to window_size
        if len(indices) < self.window_size:
            indices += [4] * (self.window_size - len(indices))
        elif len(indices) > self.window_size:
            indices = indices[:self.window_size]
        
        # Convert to tensor
        x = torch.tensor([indices], dtype=torch.long)
        
        # Encode
        with torch.no_grad():
            device = next(self.model.parameters()).device
            x = x.to(device)
            # Get embedding WITHOUT contrastive head (256D)
            # Index is 128D but model wasn't trained properly with projection
            embedding = self.model(x, use_contrastive_head=False)
            embedding = embedding.cpu().numpy()[0]
        
        return embedding
    
    def extract_seeds(self, read: str, num_seeds: int) -> List[Seed]:
        """
        Extract evenly-spaced seeds from read
        
        Args:
            read: DNA sequence
            num_seeds: Number of seeds to extract
        
        Returns:
            List of Seed objects
        """
        read_len = len(read)
        seeds = []
        
        # Calculate seed positions (evenly spaced)
        if num_seeds == 1:
            positions = [0]
        else:
            step = max(1, (read_len - self.window_size) // (num_seeds - 1))
            positions = [i * step for i in range(num_seeds)]
            # Ensure last seed doesn't exceed read length
            positions = [min(p, read_len - self.window_size) for p in positions]
        
        for pos in positions:
            # Extract window
            end = min(pos + self.window_size, read_len)
            seq = read[pos:end]
            
            # Pad if needed
            if len(seq) < self.window_size:
                seq = seq + 'N' * (self.window_size - len(seq))
            
            # Create seed
            seed = Seed(
                read_pos=pos,
                seq=seq,
                embedding=None,  # Will encode later
                candidates=[]
            )
            seeds.append(seed)
        
        return seeds
    
    def search_seed(self, seed: Seed) -> Seed:
        """
        Search seed in FAISS index and get candidates
        
        Args:
            seed: Seed to search
        
        Returns:
            Seed with candidates populated
        """
        # Encode if not already
        if seed.embedding is None:
            seed.embedding = self.encode_sequence(seed.seq)
        
        # Search in FAISS
        embedding = seed.embedding.reshape(1, -1).astype(np.float32)
        scores, indices = self.index.search(embedding, self.top_k)
        
        # Convert to candidates
        candidates = []
        for score, idx in zip(scores[0], indices[0]):
            if idx == -1:  # Invalid index
                continue
            chr_name = self.metadata['chr_names'][idx]
            pos = self.metadata['positions'][idx]
            candidates.append((chr_name, pos, float(score)))
        
        seed.candidates = candidates
        return seed
    
    def filter_seed(self, seed: Seed) -> Seed:
        """
        Classify seed as unique, ambiguous, or repeat
        
        Args:
            seed: Seed with candidates
        
        Returns:
            Seed with status updated
        """
        if len(seed.candidates) == 0:
            seed.status = "unmapped"
            return seed
        
        # Get top 2 scores
        scores = [c[2] for c in seed.candidates]
        top_score = scores[0]
        second_score = scores[1] if len(scores) > 1 else 0
        
        # Check uniqueness
        if top_score - second_score > self.uniqueness_threshold:
            seed.status = "unique"
        else:
            seed.status = "ambiguous"
        
        # Check for repeats (multiple high-scoring hits on different chromosomes)
        chr_names = [c[0] for c in seed.candidates[:10]]
        unique_chrs = len(set(chr_names))
        if unique_chrs >= 3:  # Mapping to 3+ chromosomes
            seed.status = "repeat"
        
        return seed
    
    def chain_seeds(self, seeds: List[Seed]) -> List[SeedChain]:
        """
        Chain seeds using colinearity check
        
        Args:
            seeds: List of seeds with candidates
        
        Returns:
            List of SeedChain objects sorted by score
        """
        chains = []
        
        # Group seeds by chromosome
        chr_seeds = {}
        for seed in seeds:
            for chr_name, pos, score in seed.candidates[:10]:  # Top 10 per seed
                if chr_name not in chr_seeds:
                    chr_seeds[chr_name] = []
                chr_seeds[chr_name].append((seed, pos, score))
        
        # Check each chromosome
        for chr_name, seed_list in chr_seeds.items():
            if len(seed_list) < 2:
                continue  # Need at least 2 seeds for a chain
            
            # Sort by read position
            seed_list.sort(key=lambda x: x[0].read_pos)
            
            # Check colinearity
            colinear = True
            for i in range(len(seed_list) - 1):
                seed1, pos1, score1 = seed_list[i]
                seed2, pos2, score2 = seed_list[i + 1]
                
                # Distance in read
                read_dist = seed2.read_pos - seed1.read_pos
                
                # Distance in reference
                ref_dist = pos2 - pos1
                
                # Expected distance should be similar
                error = abs(read_dist - ref_dist)
                
                if error > self.colinearity_tolerance:
                    colinear = False
                    break
            
            if colinear:
                # Create chain
                chain_seeds = [s[0] for s in seed_list]
                # FIX: Use number of anchors (like NeuralAligner) instead of sum of scores
                # This prevents bias toward longer chromosomes
                chain_score = len(seed_list)  # Count anchors (fair to all regions)
                start_pos = min(s[1] for s in seed_list)
                end_pos = max(s[1] for s in seed_list)
                
                chain = SeedChain(
                    chr_name=chr_name,
                    seeds=chain_seeds,
                    score=chain_score,
                    colinear=True,
                    start_pos=start_pos,
                    end_pos=end_pos
                )
                chains.append(chain)
        
        # Sort by score
        chains.sort(key=lambda x: x.score, reverse=True)
        
        return chains
    
    def align_read(self, read: str, read_id: str = None, return_top_k: int = 5):
        """
        Align read using adaptive seeding strategy
        
        Args:
            read: DNA sequence
            read_id: Optional read identifier
            return_top_k: Number of top candidates to return (for EXTEND phase)
                         Set to 1 for legacy mode (single best)
        
        Returns:
            List of top-k candidates (for EXTEND phase) or single dict (legacy)
        """
        # Step 1: Initial seeding (5 seeds)
        seeds = self.extract_seeds(read, num_seeds=self.min_seeds)
        
        # Step 2: Search and filter seeds
        for i, seed in enumerate(seeds):
            seeds[i] = self.search_seed(seed)
            seeds[i] = self.filter_seed(seeds[i])
        
        # Step 3: Chain seeds
        chains = self.chain_seeds(seeds)
        
        # Step 4: Decision gate with enhanced rescue logic (like NeuralAligner)
        need_rescue = False
        
        if len(chains) == 0:
            # Condition 1: No chains found
            need_rescue = True
        elif len(chains) > 0:
            # Condition 2: Best chain score too low (< 1/3 of seeds)
            # Relaxed from 1/2 for real reads with more errors
            if chains[0].score < self.min_seeds / 3:
                need_rescue = True
            
            # Condition 3: Top chains ambiguous (similar scores)
            # This indicates multi-mapping or repetitive region
            if len(chains) > 1:
                best_score = chains[0].score
                second_score = chains[1].score
                # If second-best is within 80% of best, it's ambiguous
                if second_score >= 0.8 * best_score:
                    need_rescue = True
        
        # Perform rescue seeding if needed
        if need_rescue and len(seeds) < self.max_seeds:
            # Add more seeds for better discrimination
            rescue_seeds = self.extract_seeds(read, num_seeds=self.max_seeds)
            for i, seed in enumerate(rescue_seeds):
                rescue_seeds[i] = self.search_seed(seed)
                rescue_seeds[i] = self.filter_seed(rescue_seeds[i])
            
            chains = self.chain_seeds(rescue_seeds)
        
        if len(chains) == 0:
            # Still no chains - unmapped
            return None
        
        # Return top-k candidates (sorted by seed chain score)
        top_k = min(return_top_k, len(chains))
        candidates = []
        
        for i in range(top_k):
            chain = chains[i]
            candidates.append({
                'read_id': read_id,
                'chr': chain.chr_name,
                'start': chain.start_pos,  # For EXTEND phase
                'end': chain.end_pos,      # For EXTEND phase
                'start_pos': chain.start_pos,  # Legacy compatibility
                'end_pos': chain.end_pos,      # Legacy compatibility
                'num_seeds': len(chain.seeds),
                'score': chain.score,
                'colinear': chain.colinear,
                'status': 'primary' if len(chains) == 1 else 'candidate'
            })
        
        # Legacy mode: return single best if return_top_k == 1
        if return_top_k == 1:
            return candidates[0]
        
        return candidates

if __name__ == '__main__':
    print("AdaptiveSeeder module ready")
    print("Use: seeder = AdaptiveSeeder(model, index, metadata)")
    print("     alignment = seeder.align_read(read)")
