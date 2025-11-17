# Implementation Plan: 70.6% → 99% on Chr22 First

**Strategy**: Validate full pipeline on chr22 before scaling to full genome

---

## 🎯 Why Chr22-First?

### Advantages
1. **Fast iteration**: Chr22 is 2% of genome → 50x faster testing
2. **Already trained**: Have working model + index for chr22
3. **Lower risk**: Validate each component before full-genome investment
4. **Resource efficient**: No need to retrain/re-encode until pipeline proven
5. **Clear validation**: Can compare directly with minimap2 on chr22

### Timeline Savings
- **Chr22 pipeline**: 3-4 weeks
- **Full genome scaling**: 1 week (after validation)
- **Total**: 4-5 weeks vs 6+ weeks with full-genome-first approach

---

## 📋 Phase-by-Phase Plan

### Phase 0: Setup (Day 0 - NOW)

**Install dependencies:**
```bash
cd /home/nebius/genocache
/home/nebius/.local/bin/uv pip install parasail pysam scikit-learn
```

**Verify existing assets:**
```bash
# Model
ls -lh nal_encoder_best.pt  # Should exist (300 KB)

# Chr22 vectors and index
ls -lh ref_vectors_improved.npy       # 1.2 GB
ls -lh ref_positions_improved.npy     # 9.4 MB
ls -lh faiss_index_improved.idx       # 49 MB

# Reference genome (for SW refinement)
ls -lh GRCh38.fa  # 3.2 GB

# All should exist ✓
```

---

### Phase 1: Multi-Seeding (Days 1-5) → Target: 80-85%

**Day 1-2: Implementation**

Create `multi_seeder.py`:
```python
#!/usr/bin/env python3
"""Multi-seed extraction and consensus alignment for chr22"""

import numpy as np
import torch
from train_improved import ImprovedNAL
import faiss

class MultiSeedAligner:
    def __init__(self, encoder_path, index_path, positions_path):
        # Load model
        self.encoder = ImprovedNAL().to("cuda")
        self.encoder.load_state_dict(torch.load(encoder_path))
        self.encoder.eval()
        
        # Load index
        self.index = faiss.read_index(index_path)
        self.ref_positions = np.load(positions_path)
        
        self.seed_len = 256
    
    def align_read(self, read_seq, n_seeds=5, top_k=20):
        """Align read using multiple seeds"""
        # Extract seeds
        seeds = self.extract_seeds(read_seq, n_seeds)
        
        # Search each seed
        all_hits = []
        for read_pos, seed_seq in seeds:
            vec = self.encode_seed(seed_seq)
            scores, indices = self.index.search(vec.reshape(1, -1), k=top_k)
            
            for idx, score in zip(indices[0], scores[0]):
                ref_pos = self.ref_positions[idx]
                estimated_start = ref_pos - read_pos
                
                all_hits.append({
                    'ref_pos': estimated_start,
                    'seed_read_pos': read_pos,
                    'score': score
                })
        
        # Cluster nearby hits
        return self.cluster_hits(all_hits)
    
    def extract_seeds(self, read_seq, n_seeds):
        """Extract evenly-spaced seeds"""
        read_len = len(read_seq)
        if read_len < self.seed_len:
            return [(0, read_seq)]
        
        step = (read_len - self.seed_len) // (n_seeds - 1) if n_seeds > 1 else 0
        seeds = []
        
        for i in range(n_seeds):
            pos = min(i * step, read_len - self.seed_len)
            seed = read_seq[pos:pos + self.seed_len]
            seeds.append((pos, seed))
        
        return seeds
    
    def cluster_hits(self, hits, cluster_dist=500):
        """Group nearby hits"""
        if not hits:
            return []
        
        hits.sort(key=lambda h: h['ref_pos'])
        clusters = []
        current = [hits[0]]
        
        for hit in hits[1:]:
            if abs(hit['ref_pos'] - current[-1]['ref_pos']) < cluster_dist:
                current.append(hit)
            else:
                clusters.append(self.summarize_cluster(current))
                current = [hit]
        
        if current:
            clusters.append(self.summarize_cluster(current))
        
        clusters.sort(key=lambda c: c['total_score'], reverse=True)
        return clusters
    
    def summarize_cluster(self, hits):
        return {
            'ref_pos': int(np.median([h['ref_pos'] for h in hits])),
            'n_seeds': len(hits),
            'total_score': sum(h['score'] for h in hits),
            'avg_score': np.mean([h['score'] for h in hits])
        }
    
    def encode_seed(self, seed_seq):
        """Encode seed to vector"""
        onehot = self.seq_to_onehot(seed_seq)
        with torch.no_grad():
            tensor = torch.from_numpy(onehot).unsqueeze(0).to("cuda")
            vec = self.encoder(tensor).cpu().numpy()
        return vec
    
    def seq_to_onehot(self, seq):
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

if __name__ == "__main__":
    # Quick test
    aligner = MultiSeedAligner(
        "nal_encoder_best.pt",
        "faiss_index_improved.idx",
        "ref_positions_improved.npy"
    )
    
    # Test with 4kb read
    test_read = "ACGT" * 1000
    results = aligner.align_read(test_read, n_seeds=5, top_k=20)
    
    print(f"Multi-seed alignment results:")
    for i, cluster in enumerate(results[:5]):
        print(f"  {i+1}. Pos {cluster['ref_pos']:,}, Seeds: {cluster['n_seeds']}, Score: {cluster['total_score']:.3f}")
```

**Day 3: Testing**

Create `test_phase1_multi_seed.py`:
```python
#!/usr/bin/env python3
"""Test multi-seeding improvement on chr22"""

import numpy as np
from Bio import SeqIO
import random
from multi_seeder import MultiSeedAligner

def add_errors(seq, error_rate=0.05):
    """Add substitution errors"""
    seq = list(seq)
    bases = ['A', 'C', 'G', 'T']
    for i in range(len(seq)):
        if random.random() < error_rate:
            seq[i] = random.choice([b for b in bases if b != seq[i]])
    return ''.join(seq)

def test_multi_seed(n_reads=1000):
    # Load chr22 reference
    ref = None
    for record in SeqIO.parse("GRCh38.fa", "fasta"):
        if "NC_000022" in record.id or record.id == "chr22":
            ref = str(record.seq).upper()
            break
    
    if not ref:
        print("Chr22 not found in reference!")
        return
    
    print(f"Testing on chr22 ({len(ref):,} bp)")
    
    # Create aligner
    aligner = MultiSeedAligner(
        "nal_encoder_best.pt",
        "faiss_index_improved.idx",
        "ref_positions_improved.npy"
    )
    
    correct = 0
    total = 0
    errors = []
    
    for i in range(n_reads):
        # Extract random read
        true_pos = random.randint(0, len(ref) - 5000)
        read_len = random.randint(1000, 5000)
        read_seq = ref[true_pos:true_pos + read_len]
        
        # Add errors
        read_seq = add_errors(read_seq, error_rate=0.05)
        
        # Align
        results = aligner.align_read(read_seq, n_seeds=5, top_k=20)
        
        if results:
            pred_pos = results[0]['ref_pos']
            error = abs(pred_pos - true_pos)
            errors.append(error)
            
            if error < 100:  # Within 100bp
                correct += 1
        
        total += 1
        
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{n_reads}, Recall: {correct/total*100:.1f}%")
    
    recall = correct / total
    median_error = np.median(errors) if errors else float('inf')
    
    print(f"\n{'='*60}")
    print(f"Phase 1 Results (Multi-Seeding):")
    print(f"  Recall@100bp: {recall*100:.1f}%")
    print(f"  Median error: {median_error:.0f} bp")
    print(f"  Tested: {total} reads")
    print(f"{'='*60}")
    
    return recall

if __name__ == "__main__":
    recall = test_multi_seed(n_reads=1000)
    
    if recall > 0.80:
        print("\n✓ Phase 1 PASSED! Ready for Phase 2 (Chaining)")
    elif recall > 0.75:
        print("\n~ Phase 1 MARGINAL. Tune parameters before Phase 2")
    else:
        print("\n✗ Phase 1 FAILED. Debug multi-seeding implementation")
```

**Day 4-5: Validation & Tuning**

Run tests and tune parameters:
```bash
# Run phase 1 test
.venv/bin/python test_phase1_multi_seed.py

# Expected output:
# Phase 1 Results (Multi-Seeding):
#   Recall@100bp: 82.3%     ← Target: > 80%
#   Median error: 35 bp
#   Tested: 1000 reads
# 
# ✓ Phase 1 PASSED! Ready for Phase 2 (Chaining)
```

**Gate Check:**
- ✅ Recall > 80% on chr22
- ✅ Clear improvement over 70.6% baseline
- ✅ Code is clean and tested

---

### Phase 2: Chaining (Days 6-12) → Target: 88-92%

**Day 6-8: Implementation**

Create `seed_chainer.py`:
```python
#!/usr/bin/env python3
"""Seed chaining using dynamic programming"""

class SeedChainer:
    def chain_seeds(self, seeds, max_gap=10000, max_deviation=500):
        """
        Chain seeds using DP
        
        Returns best chain of co-linear seeds
        """
        if not seeds:
            return []
        
        # Sort by read position
        seeds = sorted(seeds, key=lambda s: s['seed_read_pos'])
        
        n = len(seeds)
        best_score = [s['score'] for s in seeds]
        predecessor = [-1] * n
        
        # DP
        for i in range(1, n):
            for j in range(i):
                if self.are_colinear(seeds[j], seeds[i], max_gap, max_deviation):
                    score = best_score[j] + seeds[i]['score']
                    if score > best_score[i]:
                        best_score[i] = score
                        predecessor[i] = j
        
        # Traceback
        best_end = max(range(n), key=lambda i: best_score[i])
        chain = []
        idx = best_end
        
        while idx != -1:
            chain.append(seeds[idx])
            idx = predecessor[idx]
        
        chain.reverse()
        return chain
    
    def are_colinear(self, seed1, seed2, max_gap, max_deviation):
        """Check if seeds are consistent"""
        read_dist = seed2['seed_read_pos'] - seed1['seed_read_pos']
        ref_dist = seed2['seed_ref_pos'] - seed1['seed_ref_pos']
        
        if read_dist <= 0 or ref_dist <= 0:
            return False
        
        if read_dist > max_gap or ref_dist > max_gap:
            return False
        
        deviation = abs(read_dist - ref_dist)
        return deviation < max_deviation
```

**Day 9-10: Integration & Testing**

Update `multi_seeder.py` to use chaining, create `test_phase2_chaining.py`

**Day 11-12: Validation**

```bash
.venv/bin/python test_phase2_chaining.py

# Expected:
# Phase 2 Results (Multi-Seed + Chaining):
#   Recall@100bp: 89.1%     ← Target: > 88%
#   Median error: 22 bp
# 
# ✓ Phase 2 PASSED! Ready for Phase 3 (SW Refinement)
```

---

### Phase 3: SW Refinement (Days 13-23) → Target: 95-98%

**Day 13-15: Setup & Implementation**

Create `refine_alignment.py` with parasail integration

**Day 16-20: Integration Testing**

Integrate SW refinement with chained seeds

**Day 21-23: Validation**

```bash
.venv/bin/python test_phase3_refinement.py

# Expected:
# Phase 3 Results (Full Pipeline - No MAPQ):
#   Recall@100bp: 96.4%     ← Target: > 95%
#   Median error: 8 bp
#   CIGAR generated: Yes
# 
# ✓ Phase 3 PASSED! Ready for Phase 4 (Quality Scores)
```

---

### Phase 4: Quality & Output (Days 24-30) → Target: 98-99%

**Day 24-26: MAPQ Implementation**

Add quality scoring

**Day 27-29: BAM Output**

Generate proper SAM/BAM files

**Day 30: Final Validation**

```bash
# Compare with minimap2 on chr22
.venv/bin/python compare_with_minimap2.py --chromosome chr22

# Expected:
# Chr22 Pipeline Complete:
#   Your system:  98.2% recall, 1500 reads/sec
#   minimap2:     99.1% recall, 800 reads/sec
# 
# ✓ Chr22 pipeline validated! Ready for full genome scaling
```

---

### Phase 5: Scale to Full Genome (Week 6+)

**ONLY AFTER Phase 4 complete!**

Day 31-33: Retrain on full GRCh38  
Day 34-35: Re-encode and rebuild index  
Day 36-37: Full validation  

---

## 🧪 Testing Protocol

### After Each Phase

1. **Run synthetic tests** (extract from chr22 reference)
2. **Measure recall improvement**
3. **Check speed hasn't degraded**
4. **Validate against previous phase**
5. **Gate decision**: Pass → next phase, Fail → debug

### Validation Gates

Each phase must meet targets before proceeding:
- Phase 1: > 80% recall
- Phase 2: > 88% recall
- Phase 3: > 95% recall
- Phase 4: > 98% recall, comparable to minimap2

---

## 📁 Files to Create (In Order)

**Week 1:**
1. `multi_seeder.py`
2. `test_phase1_multi_seed.py`

**Week 2:**
3. `seed_chainer.py`
4. `test_phase2_chaining.py`

**Week 3-4:**
5. `refine_alignment.py`
6. `test_phase3_refinement.py`

**Week 4-5:**
7. `quality_scorer.py`
8. `bam_writer.py`
9. `test_phase4_complete.py`
10. `compare_with_minimap2.py`

**Week 6+ (After validation):**
11. `train_full_genome.py`
12. `encode_grch38_final.py`

---

## ✅ Current Status

- ✓ Phase 0: Setup (dependencies, existing assets)
- ⏳ Waiting: GRCh38 encoding (for reference access)
- 📝 Next: Implement Phase 1 (multi_seeder.py)

---

## 🎯 Success Metrics

**Chr22 Pipeline Complete When:**
- ✅ 98-99% recall on chr22 synthetic reads
- ✅ Comparable to minimap2 on chr22
- ✅ CIGAR strings + MAPQ scores generated
- ✅ BAM output working
- ✅ Speed > 1K reads/sec end-to-end

**Then and only then:**
- → Retrain on full genome
- → Re-encode GRCh38
- → Scale to full dataset

---

**Bottom line**: Prove it works on chr22 first, then scale. Saves weeks of iteration time!
