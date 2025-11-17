# Roadmap to 99.9% Accuracy - Complete Aligner

**Goal**: Transform neural seeder (70.6%) into complete aligner (99.9%)  
**Current**: Seeding only  
**Target**: Full alignment pipeline matching minimap2/BWA-MEM2

---

## 🎯 Architecture Overview

### Current System (Seeding Only - 70.6%)
```
Read → Neural Encoder → FAISS Search → Top-K Positions
```

### Target System (Complete Aligner - 99.9%)
```
Read → [Stage 1: Multi-Seeding]    → Seed clusters
     → [Stage 2: Chaining]         → Consistent chains
     → [Stage 3: RC Handling]      → Both strands
     → [Stage 4: Local Alignment]  → Exact alignment
     → [Stage 5: Quality Scoring]  → MAPQ + CIGAR
     → BAM Output
```

---

## 📊 Accuracy Progression Path (Chr22 Only)

| Stage | Recall | Increment | Effort | Time | Dataset |
|-------|--------|-----------|--------|------|---------|
| **Current** | 70.6% | - | - | ✅ Done | chr22 |
| + Multi-seeding | 80-85% | +10-15% | Medium | 3-5 days | chr22 |
| + RC search | 85-88% | +3-5% | Low | 1-2 days | chr22 |
| + Chaining | 88-92% | +3-5% | Medium | 5-7 days | chr22 |
| + SW/WFA refinement | 95-98% | +7-10% | High | 7-10 days | chr22 |
| + Quality scores | 98-99% | +1-2% | Medium | 3-5 days | chr22 |
| **Chr22 Complete** | **98-99%** | **+28%** | **High** | **~3-4 weeks** | chr22 |
| + Full genome retrain | 99-99.9% | +1-2% | Medium | 2-3 days | **GRCh38** |
| **Production Ready** | **99-99.9%** | **+29%** | **High** | **~4-6 weeks** | GRCh38 |

**Strategy**: Validate entire pipeline on chr22 first, then scale to full genome

---

## 🚀 Implementation Plan

### Phase 1: Multi-Seeding (Week 1) → 80-85% 

**Objective**: Extract and align multiple seeds per read

**Implementation:**
```python
class MultiSeeder:
    def extract_seeds(self, read, n_seeds=5):
        """Extract multiple seeds from different positions"""
        read_len = len(read)
        seed_positions = []
        
        # Strategy 1: Evenly spaced
        step = (read_len - SEED_LEN) // (n_seeds - 1)
        for i in range(n_seeds):
            pos = min(i * step, read_len - SEED_LEN)
            seed = read[pos:pos + SEED_LEN]
            seed_positions.append((pos, seed))
        
        return seed_positions
    
    def align_multi_seed(self, read):
        """Align using multiple seeds"""
        seeds = self.extract_seeds(read, n_seeds=5)
        all_candidates = []
        
        for read_pos, seed_seq in seeds:
            # Encode and search each seed
            vec = self.encode(seed_seq)
            scores, indices = self.index.search(vec, k=10)
            
            for rank, (idx, score) in enumerate(zip(indices, scores)):
                ref_pos = self.ref_positions[idx]
                all_candidates.append({
                    'read_pos': read_pos,
                    'ref_pos': ref_pos,
                    'score': score,
                    'rank': rank
                })
        
        # Cluster candidates by position
        return self.cluster_candidates(all_candidates)
```

**Expected improvement**: +10-15% recall  
**Effort**: 3-5 days  
**Files to create**: `multi_seeder.py`, `test_multi_seed.py`

---

### Phase 2: Reverse Complement (Week 1) → 85-88%

**Objective**: Search both forward and reverse complement strands

**Implementation:**
```python
def reverse_complement(seq):
    """Generate reverse complement"""
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]

class StrandAwareSeeder:
    def align_both_strands(self, read):
        """Search both strands"""
        # Forward strand
        fwd_seeds = self.multi_seed_align(read)
        
        # Reverse complement
        rc_read = reverse_complement(read)
        rc_seeds = self.multi_seed_align(rc_read)
        
        # Mark RC alignments
        for seed in rc_seeds:
            seed['strand'] = '-'
            seed['ref_pos_adjusted'] = self.adjust_rc_position(seed)
        
        for seed in fwd_seeds:
            seed['strand'] = '+'
        
        # Combine and rank
        all_seeds = fwd_seeds + rc_seeds
        all_seeds.sort(key=lambda x: x['score'], reverse=True)
        
        return all_seeds
```

**Expected improvement**: +3-5% recall  
**Effort**: 1-2 days  
**Critical for real ONT data** (50% of reads are RC)

---

### Phase 3: Seed Chaining (Week 2) → 88-92%

**Objective**: Link multiple seeds into consistent chains

**Implementation:**
```python
class SeedChainer:
    def chain_seeds(self, seeds, max_gap=10000, max_deviation=500):
        """
        Chain seeds that are co-linear in read and reference
        
        Good chain example:
        Read:  [seed1 @ 0bp]----[seed2 @ 500bp]----[seed3 @ 1000bp]
        Ref:   [pos1]----------[pos2]--------------[pos3]
        
        If pos2 - pos1 ≈ 500 and pos3 - pos2 ≈ 500 → Good chain!
        """
        chains = []
        
        for seed in seeds:
            # Try to extend existing chains
            extended = False
            
            for chain in chains:
                if self.is_colinear(chain[-1], seed, max_gap, max_deviation):
                    chain.append(seed)
                    extended = True
                    break
            
            if not extended:
                # Start new chain
                chains.append([seed])
        
        # Score chains
        for chain in chains:
            chain_score = self.score_chain(chain)
            chain_length = sum(s['coverage'] for s in chain)
            chain_confidence = self.chain_confidence(chain)
        
        # Return best chain
        return max(chains, key=lambda c: self.chain_score(c))
    
    def is_colinear(self, seed1, seed2, max_gap, max_deviation):
        """Check if two seeds are consistent"""
        # Distance in read
        read_dist = seed2['read_pos'] - seed1['read_pos']
        
        # Distance in reference
        ref_dist = seed2['ref_pos'] - seed1['ref_pos']
        
        # Should be similar (allowing for indels)
        deviation = abs(read_dist - ref_dist)
        
        return (0 < read_dist < max_gap and 
                0 < ref_dist < max_gap and
                deviation < max_deviation)
```

**Expected improvement**: +3-5% recall  
**Effort**: 5-7 days  
**Key for long reads**: Links evidence across read length

---

### Phase 4: Local Alignment Refinement (Week 3-4) → 95-98%

**Objective**: Precise alignment with indels

**Option A: Smith-Waterman (Simpler)**
```python
class SmithWaterman:
    def align(self, query, reference, region_start, region_end):
        """
        Standard SW alignment
        """
        # Get reference region (+/- margin for indels)
        margin = 1000
        ref_seq = self.get_reference(region_start - margin, 
                                     region_end + margin)
        
        # DP matrix
        score_matrix = self.build_matrix(query, ref_seq)
        
        # Traceback
        alignment = self.traceback(score_matrix)
        
        return {
            'cigar': alignment.cigar,
            'position': alignment.position,
            'score': alignment.score,
            'matches': alignment.matches,
            'mismatches': alignment.mismatches,
            'indels': alignment.indels
        }
```

**Option B: WFA (Faster, more complex)**
```python
# Use WFA library (Wavefront Alignment)
import pywfa

class WFAAligner:
    def __init__(self):
        self.aligner = pywfa.WavefrontAligner()
    
    def align(self, query, reference, seed_position):
        """Ultra-fast alignment using WFA"""
        # Extract reference region
        margin = 2000
        ref_region = self.get_reference(
            seed_position - margin,
            seed_position + len(query) + margin
        )
        
        # WFA alignment
        result = self.aligner.align(query, ref_region)
        
        return {
            'cigar': result.cigar,
            'position': seed_position + result.offset,
            'score': result.score
        }
```

**Expected improvement**: +7-10% recall  
**Effort**: 7-10 days  
**Critical for 99% accuracy**

---

### Phase 5: Full Genome Retraining (Week 3) → 97-99%

**Objective**: Train on diverse genomic sequences, not just chr22

**Changes:**
```python
# In train_improved.py

class FullGenomeSampler(Dataset):
    def __init__(self, ref_fa="GRCh38.fa", L=256):
        self.chromosomes = []
        
        # Load all main chromosomes
        for record in SeqIO.parse(ref_fa, "fasta"):
            if record.id.startswith('NC_'):
                self.chromosomes.append({
                    'name': record.id,
                    'seq': str(record.seq).upper(),
                    'weight': len(record.seq)  # Sample proportional to size
                })
        
        # Normalize weights
        total = sum(c['weight'] for c in self.chromosomes)
        for c in self.chromosomes:
            c['prob'] = c['weight'] / total
    
    def __getitem__(self, idx):
        # Sample chromosome proportionally
        chr_probs = [c['prob'] for c in self.chromosomes]
        chr_idx = np.random.choice(len(self.chromosomes), p=chr_probs)
        
        chrom = self.chromosomes[chr_idx]
        pos = random.randint(0, len(chrom['seq']) - self.L)
        seq = chrom['seq'][pos:pos + self.L]
        
        return seq, chr_idx, pos
```

**Expected improvement**: +2-3% recall  
**Training time**: 2-3 hours (longer than chr22)  
**When to do**: After implementing chaining

---

### Phase 6: Quality Scores & Output (Week 4-5) → 99-99.9%

**Objective**: Production-ready BAM output with MAPQ scores

**Implementation:**
```python
class QualityScorer:
    def calculate_mapq(self, alignment, alternatives):
        """
        Calculate mapping quality (Phred-scaled probability of error)
        
        MAPQ = -10 * log10(P_error)
        
        High MAPQ (>40): Unique, confident alignment
        Low MAPQ (<20): Multi-mapping, ambiguous
        """
        primary_score = alignment['score']
        
        # Check if multi-mapping
        if len(alternatives) > 1:
            second_best = alternatives[1]['score']
            score_diff = primary_score - second_best
            
            if score_diff < 0.05:  # Very close scores
                mapq = 0  # Ambiguous
            elif score_diff < 0.1:
                mapq = 10  # Low confidence
            else:
                mapq = 40  # Unique
        else:
            mapq = 60  # Single good match
        
        # Adjust by alignment identity
        identity = alignment['matches'] / alignment['length']
        mapq = int(mapq * identity)
        
        return min(60, max(0, mapq))

class BAMWriter:
    def write_alignment(self, read, alignment):
        """Write SAM/BAM record"""
        import pysam
        
        sam_record = pysam.AlignedSegment()
        sam_record.query_name = read.id
        sam_record.query_sequence = str(read.seq)
        sam_record.reference_id = alignment['chrom']
        sam_record.reference_start = alignment['position']
        sam_record.mapping_quality = alignment['mapq']
        sam_record.cigar = alignment['cigar']
        sam_record.flag = 16 if alignment['strand'] == '-' else 0
        
        return sam_record
```

**Expected improvement**: +1-2% (edge cases, quality filtering)  
**Effort**: 3-5 days  
**Required for production use**

---

## 🔧 Detailed Implementation Steps

### Step 1: Install Additional Dependencies (Now)

```bash
cd /home/nebius/genocache

# Smith-Waterman implementation
/home/nebius/.local/bin/uv pip install parasail

# Or WFA (faster, harder to integrate)
# git clone https://github.com/smarco/WFA2-lib
# cd WFA2-lib && make && cd bindings/python && pip install .

# BAM/SAM handling
/home/nebius/.local/bin/uv pip install pysam

# For testing
/home/nebius/.local/bin/uv pip install scikit-learn
```

### Step 2: Implement Multi-Seeder (Priority 1)

**Create**: `multi_seeder.py`

```python
#!/usr/bin/env python3
"""
Multi-seed extraction and consensus alignment
"""
import numpy as np
import torch
from collections import defaultdict

class MultiSeedAligner:
    def __init__(self, encoder, index, ref_positions, ref_chromosomes):
        self.encoder = encoder
        self.index = index
        self.ref_positions = ref_positions
        self.ref_chromosomes = ref_chromosomes
        self.seed_len = 256
    
    def extract_seeds(self, read_seq, n_seeds=5, strategy='evenly_spaced'):
        """
        Extract multiple seeds from a read
        
        Strategies:
        - 'evenly_spaced': Seeds at regular intervals
        - 'high_quality': Seeds from high-quality regions (low N content)
        - 'random': Random positions
        """
        read_len = len(read_seq)
        
        if read_len < self.seed_len:
            return [(0, read_seq)]
        
        seeds = []
        
        if strategy == 'evenly_spaced':
            # Evenly distribute seeds across read
            if read_len < self.seed_len * n_seeds:
                # Overlapping seeds
                step = (read_len - self.seed_len) // (n_seeds - 1) if n_seeds > 1 else 0
            else:
                # Non-overlapping seeds
                step = (read_len - self.seed_len) // n_seeds
            
            for i in range(n_seeds):
                pos = i * step
                if pos + self.seed_len <= read_len:
                    seed = read_seq[pos:pos + self.seed_len]
                    seeds.append((pos, seed))
        
        elif strategy == 'high_quality':
            # Prefer regions with fewer Ns
            for i in range(0, read_len - self.seed_len + 1, self.seed_len // 2):
                seed = read_seq[i:i + self.seed_len]
                n_count = seed.count('N')
                quality = (self.seed_len - n_count) / self.seed_len
                seeds.append((i, seed, quality))
            
            # Take top N by quality
            seeds.sort(key=lambda x: x[2], reverse=True)
            seeds = [(pos, seq) for pos, seq, qual in seeds[:n_seeds]]
        
        return seeds
    
    def align_read(self, read_seq, top_k=20, min_score=0.5):
        """
        Align read using multiple seeds
        
        Returns clustered candidate positions
        """
        # Extract seeds
        seeds = self.extract_seeds(read_seq, n_seeds=5)
        
        # Align each seed
        all_hits = []
        
        for read_pos, seed_seq in seeds:
            # Encode seed
            vec = self.encode_seed(seed_seq)
            
            # Search index
            scores, indices = self.index.search(vec.reshape(1, -1), k=top_k)
            
            # Store hits
            for rank, (idx, score) in enumerate(zip(indices[0], scores[0])):
                if score < min_score:
                    continue
                
                ref_pos = self.ref_positions[idx]
                chrom = self.ref_chromosomes[idx]
                
                # Adjust reference position based on seed position in read
                estimated_read_start = ref_pos - read_pos
                
                all_hits.append({
                    'chrom': chrom,
                    'ref_pos': estimated_read_start,
                    'seed_read_pos': read_pos,
                    'seed_ref_pos': ref_pos,
                    'score': score,
                    'rank': rank
                })
        
        # Cluster hits by position
        return self.cluster_hits(all_hits, cluster_dist=500)
    
    def cluster_hits(self, hits, cluster_dist=500):
        """
        Group hits that are nearby (likely same alignment)
        """
        if not hits:
            return []
        
        # Sort by chromosome and position
        hits.sort(key=lambda h: (h['chrom'], h['ref_pos']))
        
        clusters = []
        current_cluster = [hits[0]]
        
        for hit in hits[1:]:
            last = current_cluster[-1]
            
            # Same chromosome and nearby?
            if (hit['chrom'] == last['chrom'] and 
                abs(hit['ref_pos'] - last['ref_pos']) < cluster_dist):
                current_cluster.append(hit)
            else:
                # Save cluster and start new one
                clusters.append(self.summarize_cluster(current_cluster))
                current_cluster = [hit]
        
        # Don't forget last cluster
        if current_cluster:
            clusters.append(self.summarize_cluster(current_cluster))
        
        # Sort by total score
        clusters.sort(key=lambda c: c['total_score'], reverse=True)
        return clusters
    
    def summarize_cluster(self, hits):
        """Summarize a cluster of hits"""
        return {
            'chrom': hits[0]['chrom'],
            'ref_pos': int(np.median([h['ref_pos'] for h in hits])),
            'n_seeds': len(hits),
            'total_score': sum(h['score'] for h in hits),
            'avg_score': np.mean([h['score'] for h in hits]),
            'hits': hits
        }
    
    def encode_seed(self, seed_seq):
        """Encode a seed sequence to vector"""
        # Convert to one-hot
        onehot = self.seq_to_onehot(seed_seq)
        
        # Encode with neural network
        with torch.no_grad():
            tensor = torch.from_numpy(onehot).unsqueeze(0).to(next(self.encoder.parameters()).device)
            vec = self.encoder(tensor).cpu().numpy()
        
        return vec
    
    def seq_to_onehot(self, seq):
        """Convert sequence to one-hot encoding"""
        import random
        arr = np.zeros((len(seq), 4), dtype=np.float32)
        for i, ch in enumerate(seq.upper()):
            if ch == "A": arr[i,0] = 1
            elif ch == "C": arr[i,1] = 1
            elif ch == "G": arr[i,2] = 1
            elif ch == "T": arr[i,3] = 1
            elif ch == "N": arr[i, random.randrange(4)] = 0.25
            else: arr[i, random.randrange(4)] = 1.0
        return arr
```

**Expected**: Multi-seeding alone will improve 70.6% → 80-85%

---

### Step 3: Test Multi-Seeder (Priority 1)

**Create**: `test_multi_seed.py`

```python
#!/usr/bin/env python3
"""Test multi-seeding improvement"""

import numpy as np
import torch
from multi_seeder import MultiSeedAligner
from train_improved import ImprovedNAL
import faiss

# Load model and index
model = ImprovedNAL().to("cuda")
model.load_state_dict(torch.load("nal_encoder_best.pt"))
model.eval()

index = faiss.read_index("faiss_index_improved.idx")
ref_positions = np.load("ref_positions_improved.npy")
ref_chromosomes = np.zeros(len(ref_positions), dtype=np.int32)  # Chr22 = 0

# Create multi-seeder
aligner = MultiSeedAligner(model, index, ref_positions, ref_chromosomes)

# Test read
test_read = "ACGT" * 1000  # 4kb read

# Align with multi-seeding
results = aligner.align_read(test_read, top_k=20)

print(f"Found {len(results)} candidate alignments:")
for i, cluster in enumerate(results[:5]):
    print(f"  {i+1}. Chr {cluster['chrom']}, pos {cluster['ref_pos']:,}")
    print(f"     Seeds: {cluster['n_seeds']}, Score: {cluster['total_score']:.3f}")
```

---

### Step 4: Implement Chaining (Priority 2)

**Create**: `seed_chainer.py`

```python
#!/usr/bin/env python3
"""Seed chaining for co-linear alignment"""

class SeedChainer:
    def chain_seeds(self, seeds, max_gap=10000, max_deviation=500):
        """
        Chain seeds using dynamic programming
        
        Similar to minimap2's chaining algorithm
        """
        if not seeds:
            return []
        
        # Sort seeds by read position
        seeds = sorted(seeds, key=lambda s: s['seed_read_pos'])
        
        n = len(seeds)
        # DP: best_score[i] = best chain ending at seed i
        best_score = [s['score'] for s in seeds]
        predecessor = [-1] * n
        
        # DP: find best chain
        for i in range(1, n):
            for j in range(i):
                # Check if seeds i and j are co-linear
                if self.are_colinear(seeds[j], seeds[i], max_gap, max_deviation):
                    score = best_score[j] + seeds[i]['score']
                    if score > best_score[i]:
                        best_score[i] = score
                        predecessor[i] = j
        
        # Traceback to get best chain
        best_end = max(range(n), key=lambda i: best_score[i])
        chain = []
        idx = best_end
        
        while idx != -1:
            chain.append(seeds[idx])
            idx = predecessor[idx]
        
        chain.reverse()
        return chain
    
    def are_colinear(self, seed1, seed2, max_gap, max_deviation):
        """Check if two seeds are consistent"""
        # Must be same chromosome
        if seed1.get('chrom') != seed2.get('chrom'):
            return False
        
        # Distance in read
        read_dist = seed2['seed_read_pos'] - seed1['seed_read_pos']
        
        # Distance in reference
        ref_dist = seed2['seed_ref_pos'] - seed1['seed_ref_pos']
        
        # Both should be positive and similar
        if read_dist <= 0 or ref_dist <= 0:
            return False
        
        if read_dist > max_gap or ref_dist > max_gap:
            return False
        
        # Allow some deviation (indels)
        deviation = abs(read_dist - ref_dist)
        return deviation < max_deviation
```

**Expected**: Chaining will improve 80% → 88-92%

---

### Step 5: Integrate Smith-Waterman (Priority 3)

**Create**: `refine_alignment.py`

```python
#!/usr/bin/env python3
"""Local alignment refinement with Smith-Waterman"""

import parasail
from Bio import SeqIO

class AlignmentRefiner:
    def __init__(self, ref_fa="GRCh38.fa"):
        # Load reference (memory-mapped for efficiency)
        self.reference = {}
        for record in SeqIO.parse(ref_fa, "fasta"):
            self.reference[record.id] = str(record.seq).upper()
    
    def refine_alignment(self, read_seq, chrom, seed_pos, margin=2000):
        """
        Refine alignment using Smith-Waterman
        
        Args:
            read_seq: Query sequence
            chrom: Chromosome ID
            seed_pos: Approximate position from seeding
            margin: Bases to include around seed position
        
        Returns:
            Refined alignment with CIGAR string
        """
        # Extract reference region
        ref_seq = self.reference[chrom]
        ref_start = max(0, seed_pos - margin)
        ref_end = min(len(ref_seq), seed_pos + len(read_seq) + margin)
        ref_region = ref_seq[ref_start:ref_end]
        
        # Smith-Waterman alignment
        result = parasail.sw_trace_striped_16(
            read_seq, 
            ref_region,
            open=5,      # Gap open penalty
            extend=2,    # Gap extension penalty
            matrix=parasail.blosum62  # Scoring matrix
        )
        
        # Parse CIGAR
        cigar = result.cigar.decode
        
        # Calculate alignment position
        aligned_pos = ref_start + result.ref_begin
        
        return {
            'position': aligned_pos,
            'cigar': self.format_cigar(result.cigar),
            'score': result.score,
            'matches': result.matches,
            'length': result.len_query,
            'identity': result.matches / result.len_query if result.len_query > 0 else 0
        }
    
    def format_cigar(self, cigar_obj):
        """Convert parasail cigar to SAM format"""
        cigar_str = ""
        for op, length in cigar_obj:
            cigar_str += f"{length}{op}"
        return cigar_str
```

**Expected**: SW refinement will improve 88% → 95-98%

---

## 📅 Implementation Timeline (Chr22-First Strategy)

### Week 1: Multi-Seeding + RC on Chr22 (80-88%)
- Day 1-2: Implement `multi_seeder.py`
- Day 3: Test on chr22, validate improvement
- Day 4: Add reverse complement handling
- Day 5: Integration testing **on chr22 only**
- **Validation**: 70.6% → 80-85% on chr22 synthetic reads

### Week 2: Chaining on Chr22 (88-92%)
- Day 1-3: Implement `seed_chainer.py`
- Day 4: Dynamic programming optimization
- Day 5: Test and validate **on chr22 only**
- **Validation**: 80% → 88-92% on chr22 synthetic reads

### Week 3-4: Local Alignment on Chr22 (95-98%)
- Day 1-2: Install and test parasail/WFA
- Day 3-5: Integrate SW refinement
- Day 6-7: Test on full reads **from chr22**
- Day 8-10: Debug and optimize **using chr22**
- **Validation**: 88% → 95-98% on chr22 synthetic + real reads

### Week 4-5: Quality & Output on Chr22 (98-99%)
- Day 1-2: Implement MAPQ calculation
- Day 3-4: BAM/SAM output
- Day 5: End-to-end testing **on chr22**
- Day 6-7: Bug fixes and optimization
- **Validation**: 95% → 98-99% on chr22, compare with minimap2

### Week 6+: Scale to Full Genome (ONLY AFTER CHR22 VALIDATED)
- Day 1: Prepare full genome training data
- Day 2-3: Retrain model on full GRCh38
- Day 4-5: Re-encode GRCh38, rebuild full index
- Day 6-7: Test on all chromosomes
- **Final Target**: 99-99.9% on full genome

---

## 🧪 Testing Strategy

### After Each Phase

**Test Suite**: `test_pipeline.py`

```python
def test_accuracy(aligner, n_reads=1000):
    """Test accuracy at each stage"""
    from Bio import SeqIO
    import random
    
    # Extract test reads from reference
    ref = SeqIO.read("GRCh38.fa", "fasta")
    
    correct = 0
    total = 0
    
    for i in range(n_reads):
        # Extract random read
        true_pos = random.randint(0, len(ref) - 5000)
        read_len = random.randint(1000, 10000)
        read_seq = str(ref.seq[true_pos:true_pos + read_len])
        
        # Add sequencing errors (5%)
        read_seq = add_errors(read_seq, error_rate=0.05)
        
        # Align
        result = aligner.align(read_seq)
        
        # Check if correct
        if abs(result['position'] - true_pos) < 100:
            correct += 1
        
        total += 1
    
    recall = correct / total
    print(f"Recall: {recall*100:.1f}%")
    return recall
```

Run after each phase to track improvement.

---

## 🎯 Success Criteria (Chr22-First)

### Phase 1 Complete (Multi-Seeding on Chr22)
- ✅ Recall@1 > 80% on chr22
- ✅ Median error < 50 bp on chr22
- ✅ Speed > 10K reads/sec
- **Gate**: Must validate before proceeding to Phase 2

### Phase 2 Complete (RC + Chaining on Chr22)
- ✅ Recall@1 > 88% on chr22
- ✅ Handles both strands correctly
- ✅ Chains seeds co-linearly
- **Gate**: Must validate before proceeding to Phase 3

### Phase 3 Complete (SW Refinement on Chr22)
- ✅ Recall@1 > 95% on chr22
- ✅ Median error < 10 bp on chr22
- ✅ CIGAR strings generated
- **Gate**: Must validate before proceeding to Phase 4

### Phase 4 Complete (Full Chr22 Pipeline)
- ✅ Recall@1 > 98% on chr22
- ✅ MAPQ scores calculated
- ✅ BAM output generated
- ✅ Speed > 1K reads/sec (end-to-end)
- ✅ Matches minimap2 accuracy on chr22
- **Gate**: Must validate before scaling to full genome

### Phase 5 Complete (Full Genome - FINAL)
- ✅ Recall@1 > 99% on all chromosomes
- ✅ Speed > 1K reads/sec on full genome
- ✅ Production-ready BAM output
- ✅ Comparable to minimap2 on HG002 full dataset

---

## 📊 Benchmarking Plan

**Compare against minimap2:**

```bash
# Your system
time .venv/bin/python align_complete.py \
  --reads hg002_10k.fastq \
  --ref GRCh38.fa \
  --out genocache.bam

# minimap2
time minimap2 -ax map-ont GRCh38.fa hg002_10k.fastq > minimap2.bam

# Compare accuracy
.venv/bin/python compare_aligners.py \
  --truth minimap2.bam \
  --test genocache.bam
```

**Metrics to track:**
- Recall (% reads aligned)
- Precision (% correct alignments)
- Median error (bp from true position)
- Speed (reads/sec)
- MAPQ correlation with minimap2

---

## 🚀 Quick Start - Chr22 Validation Flow

### Phase 0: Setup (Now)

```bash
# 1. Install dependencies
/home/nebius/.local/bin/uv pip install parasail pysam scikit-learn

# 2. Wait for GRCh38 encoding to complete (for reference genome access)
# 3. Use existing chr22 model and index (already have 70.6% baseline)
```

### Phase 1: Multi-Seeding on Chr22 (Week 1) → Target: 80%

```bash
# 1. Create multi_seeder.py (already provided in roadmap)
# 2. Test on chr22 synthetic reads
.venv/bin/python test_multi_seed.py

# 3. Validate improvement
.venv/bin/python test_improved_alignment.py --mode multi_seed

# 4. Expected: 70.6% → 80-85% on chr22
# 5. If validated → proceed to Phase 2
```

### Phase 2: Add Chaining on Chr22 (Week 2) → Target: 90%

```bash
# 1. Implement seed_chainer.py
# 2. Test combined multi-seed + chaining on chr22
.venv/bin/python test_chaining.py

# 3. Expected: 80% → 88-92% on chr22
# 4. If validated → proceed to Phase 3
```

### Phase 3: Add SW Refinement on Chr22 (Week 3-4) → Target: 95-99%

```bash
# 1. Integrate Smith-Waterman refinement
# 2. End-to-end testing on chr22
.venv/bin/python test_pipeline.py --chromosome chr22

# 3. Expected: 88% → 95-98% on chr22
# 4. If validated → proceed to Phase 4
```

### Phase 4: Full Chr22 Pipeline (Week 4-5) → Target: 98-99%

```bash
# 1. Add quality scores and BAM output
# 2. Complete end-to-end pipeline on chr22
# 3. Compare against minimap2 on chr22 reads
.venv/bin/python compare_aligners.py --region chr22

# 4. Expected: 98-99% on chr22
# 5. If validated → THEN scale to full genome
```

### Phase 5: Scale to Full Genome (Week 6+) → Target: 99-99.9%

**Only after chr22 pipeline is validated!**

```bash
# 1. Retrain model on full GRCh38 (not just chr22)
.venv/bin/python train_full_genome.py

# 2. Re-encode full genome (we'll have infrastructure ready)
.venv/bin/python encode_grch38_robust.py

# 3. Test on all chromosomes
.venv/bin/python test_pipeline.py --all-chromosomes

# 4. Final validation against minimap2
```

---

## 📁 Files to Create

```
genocache/
├── multi_seeder.py          ← Phase 1
├── test_multi_seed.py       ← Phase 1
├── reverse_complement.py    ← Phase 1
├── seed_chainer.py          ← Phase 2
├── test_chaining.py         ← Phase 2
├── refine_alignment.py      ← Phase 3
├── quality_scorer.py        ← Phase 4
├── bam_writer.py            ← Phase 4
├── align_complete.py        ← Full pipeline
├── test_pipeline.py         ← Testing
├── compare_aligners.py      ← Benchmarking
└── train_full_genome.py     ← Retraining
```

---

## 💡 Key Insights

### Why 70% → 99% is Achievable

**1. Neural seeding is already working**
- 70% proves the embedding space is meaningful
- FAISS search is fast and accurate
- Foundation is solid

**2. Traditional methods are well-understood**
- Chaining: Standard DP algorithm
- SW/WFA: Mature libraries available
- MAPQ: Established formulas

**3. Incremental approach de-risks**
- Each phase adds 5-15% accuracy
- Can test and validate at each step
- Early phases are simpler

**4. Timeline is realistic**
- 4-6 weeks for full implementation
- Each component is 1-2 weeks of work
- Plenty of examples to reference (minimap2, BWA-MEM)

---

## ⚠️ Potential Challenges

### Challenge 1: Speed vs Accuracy Tradeoff
**Issue**: Smith-Waterman is slow  
**Solution**: Use WFA (10-100x faster), or only refine top candidates

### Challenge 2: Parameter Tuning
**Issue**: Many parameters (gap penalties, chain thresholds, etc.)  
**Solution**: Use minimap2 defaults as starting point

### Challenge 3: Edge Cases
**Issue**: Chimeric reads, low quality regions, repeats  
**Solution**: Add quality filters, handle separately

### Challenge 4: Memory for Full Genome
**Issue**: GRCh38 index is 40-50 GB  
**Solution**: Already using FAISS compression, works on chr22

---

## 🎓 Learning Resources

**Chaining algorithms:**
- minimap2 paper: https://arxiv.org/abs/1708.01492
- "Fast and accurate genomic analyses using genome graphs"

**Smith-Waterman:**
- Parasail library: https://github.com/jeffdaily/parasail
- WFA paper: https://doi.org/10.1093/bioinformatics/btaa777

**MAPQ calculation:**
- BWA paper: https://arxiv.org/abs/1303.3997
- minimap2 source code (mapq calculation)
