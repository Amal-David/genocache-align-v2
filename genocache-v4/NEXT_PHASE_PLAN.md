# GenoCache V4 - Next Phase: Adaptive Chaining + WFA GPU

**Status:** Waiting for current benchmark to complete  
**Start:** After index building + minimap2 benchmark finishes  
**Goal:** Implement full NeuralAligner-style pipeline with GPU-accelerated alignment

---

## Phase Overview

Current state: We have neural seeding working (96.6% accuracy)  
Next step: Add adaptive chaining + precise alignment for complete pipeline

**Two Major Components:**

1. **Adaptive Chaining Strategy** (from NeuralAligner)
   - Start with 5 seeds per read
   - Add rescue seeds (up to 16) if ambiguous
   - Chain seeds using dynamic programming
   - Filter repeats and low-confidence regions

2. **WFA GPU** (Wavefront Alignment)
   - Replace slow CPU alignment with GPU WFA
   - Precise base-level alignment after seeding
   - Handle indels and complex variants
   - Real data validation (actual ONT reads, not synthetic)

---

## 1. Adaptive Chaining Strategy (NeuralAligner Implementation)

### 1.1 Overview

**Current approach (simple):**
- Extract 1 seed (first 512bp)
- Search index for top match
- Return best position

**NeuralAligner approach (adaptive):**
- Extract 5 seeds (spaced evenly across read)
- Search each seed independently
- Chain seeds using DP (dynamic programming)
- Add rescue seeds if chaining fails
- Filter ambiguous/repetitive regions

### 1.2 Detailed Algorithm

```
Input: DNA read (length L, e.g. 10kb ONT read)

STEP 1: Initial Seeding (5 seeds)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Read:    [===================== 10kb read ====================]
Seeds:    ^          ^          ^          ^          ^
         pos 0     pos 2.5k   pos 5k    pos 7.5k   pos 10k-512

For each seed position i:
  - Extract 512bp window at position i
  - Encode to 128D embedding
  - FAISS search for top-k=32 matches
  - Store candidates with scores

Result: 5 × 32 = 160 candidate positions


STEP 2: Seed Filtering
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

For each seed i:
  candidates = top-32 matches
  
  # Check uniqueness
  if max_score - second_score > threshold:
    seed[i].status = "unique"
    seed[i].best = max_score position
  else:
    seed[i].status = "ambiguous"
    seed[i].candidates = top-10 positions
    
  # Check for repeats
  if many high-scoring matches on different chromosomes:
    seed[i].status = "repeat"


STEP 3: Seed Chaining (Dynamic Programming)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Goal: Find consistent path through seeds

For each chromosome c:
  # Collect seeds mapping to chromosome c
  seeds_on_c = [seed for seed in seeds if seed.chr == c]
  
  if len(seeds_on_c) < 2:
    continue  # Not enough evidence
  
  # Sort by read position
  seeds_on_c.sort(by read_pos)
  
  # Check colinearity (read order = genome order)
  colinear = True
  for i in range(len(seeds_on_c) - 1):
    read_dist = seeds_on_c[i+1].read_pos - seeds_on_c[i].read_pos
    ref_dist = seeds_on_c[i+1].ref_pos - seeds_on_c[i].ref_pos
    
    expected_dist = ref_dist
    error = abs(read_dist - expected_dist)
    
    if error > tolerance:  # e.g., 1kb
      colinear = False
      break
  
  if colinear:
    # Compute chain score
    score = sum(seed.score for seed in seeds_on_c)
    chains.append((c, seeds_on_c, score))

# Select best chain
best_chain = max(chains, key=lambda x: x[2])


STEP 4: Decision Gate
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Case 1: Strong unique chain found (5/5 seeds agree)
  → Accept alignment
  → Proceed to WFA for precise alignment

Case 2: Weak chain (2-4 seeds agree)
  → Add rescue seeds (up to 11 more)
  → Re-run chaining
  → If still weak, mark as "low quality"

Case 3: No chain found (all seeds disagree)
  → Try dense seeding (16 seeds)
  → If still fails, mark as "unmapped"

Case 4: Multiple equally good chains
  → Mark as "ambiguous" (multi-mapper)
  → Report secondary alignments


STEP 5: Gap Filling (Optional)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

If large gaps between seeds:
  - Add intermediate seeds in gap regions
  - Useful for detecting structural variants
```

### 1.3 Implementation Files

**File 1: `adaptive_seeding.py`**
```python
class AdaptiveSeeder:
    def __init__(self, model, index, min_seeds=5, max_seeds=16):
        self.model = model
        self.index = index
        self.min_seeds = min_seeds
        self.max_seeds = max_seeds
    
    def extract_seeds(self, read, num_seeds):
        """Extract evenly-spaced seeds from read"""
        
    def search_seeds(self, seeds):
        """Search each seed in FAISS index"""
        
    def filter_seeds(self, seed_results):
        """Classify seeds as unique/ambiguous/repeat"""
        
    def chain_seeds(self, seed_results):
        """DP-based chaining"""
        
    def decide_rescue(self, chain_result):
        """Determine if rescue seeds needed"""
        
    def align_read(self, read):
        """Main entry point: adaptive seeding pipeline"""
```

**File 2: `seed_chaining.py`**
```python
class SeedChainer:
    def __init__(self, gap_tolerance=1000, min_chain_length=2):
        self.gap_tolerance = gap_tolerance
        self.min_chain_length = min_chain_length
    
    def check_colinearity(self, seeds):
        """Check if seeds form colinear chain"""
        
    def compute_chain_score(self, seeds):
        """Score chain quality"""
        
    def dynamic_programming_chain(self, seeds):
        """DP algorithm for optimal chaining"""
        
    def filter_chains(self, chains):
        """Filter overlapping/weak chains"""
```

### 1.4 Expected Improvements

**Before (current):**
- 1 seed per read
- Simple top-match selection
- No validation of correctness
- Accuracy: 96.6%

**After (adaptive):**
- 5-16 seeds per read
- Chaining validation
- Rescue seeding for hard cases
- Expected accuracy: **98-99%+** (NeuralAligner level)

---

## 2. WFA GPU (Wavefront Alignment)

### 2.1 Overview

**What is WFA?**
- Wavefront Alignment Algorithm (BiWFA, 2020)
- Exact alignment (not heuristic like Smith-Waterman)
- O(s²) time where s = edit distance (not O(nm) like SW)
- Very fast for high-identity sequences (ONT: 90-95%)

**Why GPU?**
- WFA is highly parallelizable
- Batch many alignments simultaneously
- 10-100× speedup vs CPU

### 2.2 WFA Algorithm Basics

```
Input: Query sequence (read), Reference sequence (genome region)

WFA uses "wavefronts" instead of full DP matrix:
  - Wavefront k: all positions at edit distance k
  - Expand wavefronts iteratively until reaching end
  - Much faster than filling entire DP matrix

Example:
  Query:    ACGTACGT
  Reference: ACGAACGT
              |||*||||
  
  Wavefront 0: Perfect matches (first 3 bp)
  Wavefront 1: After 1 error (skip 'T' in query or 'A' in ref)
  Continue until end reached
  
  Result: Edit distance = 1, alignment path
```

### 2.3 GPU Implementation Options

**Option 1: Use existing WFA-GPU library**
- Repo: https://github.com/quim0/WFA-GPU
- Pros: Already implemented, tested
- Cons: May need integration work

**Option 2: Use NVIDIA Clara Parabricks**
- Commercial GPU alignment toolkit
- Pros: Highly optimized
- Cons: Proprietary, licensing

**Option 3: Custom CUDA kernel**
- Write our own WFA GPU implementation
- Pros: Full control, optimization
- Cons: Time-consuming, complex

**Recommendation: Option 1 (WFA-GPU library)**
- Open source, well-maintained
- Proven performance
- Can customize if needed

### 2.4 Integration Plan

```
Current pipeline:
  Read → Neural seed → Top match → Position

New pipeline:
  Read → Neural seeds (5-16) → Chaining → Candidate region → WFA GPU → Precise alignment

Detailed:
  1. Neural seeding identifies region (e.g., chr1:1000000-1010000)
  2. Extract reference sequence from that region
  3. Batch multiple alignments together
  4. Run WFA GPU for precise base-level alignment
  5. Output CIGAR string (alignment details)
```

### 2.5 Implementation Files

**File 1: `wfa_gpu_wrapper.py`**
```python
import wfagpu  # WFA-GPU library

class WFAGPUAligner:
    def __init__(self, batch_size=256):
        self.batch_size = batch_size
        self.wfa = wfagpu.WFAAligner()
    
    def align_batch(self, queries, references):
        """Batch alignment using WFA GPU"""
        
    def parse_cigar(self, alignment):
        """Convert WFA output to CIGAR string"""
        
    def compute_identity(self, cigar):
        """Calculate % identity from CIGAR"""
```

**File 2: `complete_pipeline.py`**
```python
class CompletePipeline:
    def __init__(self, model, index, genome, wfa_aligner):
        self.seeder = AdaptiveSeeder(model, index)
        self.chainer = SeedChainer()
        self.wfa = wfa_aligner
        self.genome = genome
    
    def align_read(self, read):
        # Step 1: Adaptive seeding
        seeds = self.seeder.extract_seeds(read, num=5)
        seed_results = self.seeder.search_seeds(seeds)
        
        # Step 2: Chaining
        chains = self.chainer.chain_seeds(seed_results)
        
        if not chains:
            # Rescue seeding
            seeds = self.seeder.extract_seeds(read, num=16)
            seed_results = self.seeder.search_seeds(seeds)
            chains = self.chainer.chain_seeds(seed_results)
        
        if not chains:
            return None  # Unmapped
        
        best_chain = chains[0]
        
        # Step 3: Extract reference region
        chr_name = best_chain.chr
        start = best_chain.start - 1000  # Buffer
        end = best_chain.end + 1000
        ref_seq = self.genome.fetch(chr_name, start, end)
        
        # Step 4: WFA GPU alignment
        alignment = self.wfa.align(read, ref_seq)
        
        # Step 5: Format output
        return {
            'chr': chr_name,
            'pos': start + alignment.ref_start,
            'cigar': alignment.cigar,
            'mapq': alignment.quality,
            'identity': alignment.identity
        }
```

### 2.6 Expected Performance

**Current (seeding only):**
- Output: Approximate position (±1kb)
- Time: 1.5 ms/read
- Throughput: 700 reads/sec

**After WFA GPU:**
- Output: Precise alignment with CIGAR
- Time: 1.5 ms (seed) + 0.5 ms (WFA) = 2 ms/read
- Throughput: 500 reads/sec
- Accuracy: Base-level (not just region)

**Comparison:**
- minimap2: 5-10 ms/read, 100-200 reads/sec
- GenoCache + WFA GPU: 2 ms/read, 500 reads/sec
- **Speedup: 2.5-5× faster!**

---

## 3. Real Data Validation

### 3.1 Test Datasets

**Dataset 1: ONT chr22 (Real Data)**
- Source: Human HG002 chr22 reads
- Error rate: ~5-10% (typical ONT)
- Read length: 10-50kb
- Quantity: 10,000 reads

**Dataset 2: ONT whole genome**
- Source: Human HG002 whole genome
- Full coverage test
- Quantity: 100,000 reads

**Dataset 3: PacBio HiFi (High accuracy)**
- Source: Human HG002 HiFi
- Error rate: ~1% (very high quality)
- Read length: 10-20kb
- Test model generalization

### 3.2 Validation Metrics

**Metric 1: Alignment Accuracy**
- Compare vs minimap2 alignments
- Check agreement on chr, position, strand
- Tolerance: ±100bp (after WFA)

**Metric 2: CIGAR String Quality**
- Compare CIGAR strings vs minimap2
- Check edit distance agreement
- Validate indel detection

**Metric 3: Identity Distribution**
- Histogram of % identity
- Should match ONT characteristics (85-95%)

**Metric 4: Speed Comparison**
- GenoCache + WFA GPU vs minimap2
- Measure end-to-end time
- Report reads/second

**Metric 5: Variant Detection**
- Can we detect SNPs?
- Can we detect indels?
- Can we detect structural variants?

### 3.3 Real Data Pipeline

```bash
# Step 1: Download real ONT data
wget <HG002 ONT reads URL>

# Step 2: Run GenoCache pipeline
python complete_pipeline.py \
  --model fullgenome_best_sep12.4035_epoch30.pt \
  --index genocache_v4_production.index \
  --genome GRCh38.fa \
  --reads hg002_ont.fastq \
  --output genocache_alignments.sam

# Step 3: Run minimap2 (ground truth)
minimap2 -ax map-ont GRCh38.fa hg002_ont.fastq > minimap2_alignments.sam

# Step 4: Compare results
python compare_alignments.py \
  --geocache genocache_alignments.sam \
  --minimap2 minimap2_alignments.sam \
  --output comparison_report.json
```

---

## 4. Implementation Timeline

**Phase 1: Adaptive Chaining (2-3 days)**
- Day 1: Implement adaptive seeding
- Day 2: Implement seed chaining
- Day 3: Test and validate on synthetic data

**Phase 2: WFA GPU Integration (2-3 days)**
- Day 1: Setup WFA-GPU library
- Day 2: Integrate with pipeline
- Day 3: Test and optimize batching

**Phase 3: Real Data Validation (1-2 days)**
- Day 1: Run on HG002 ONT data
- Day 2: Analyze results, compare vs minimap2

**Total: ~1 week**

---

## 5. Success Criteria

**Minimum Success:**
- Adaptive chaining increases accuracy from 96.6% → 98%
- WFA GPU integration works
- Real data validation shows comparable accuracy to minimap2
- Speed: 2-3× faster than minimap2

**Target Success:**
- Accuracy: 99%+ (NeuralAligner level)
- Speed: 4-5× faster than minimap2
- Works on real ONT data (not just synthetic)
- Complete SAM/BAM output

**Stretch Goals:**
- Better than minimap2 accuracy in some cases
- 10× speedup with optimizations
- Support for different read types (PacBio, Illumina)
- Structural variant detection

---

## 6. Dependencies & Setup

**Software Dependencies:**
```bash
# WFA-GPU library
git clone https://github.com/quim0/WFA-GPU
cd WFA-GPU
mkdir build && cd build
cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
make -j8

# Python bindings
pip install pywfa

# SAM/BAM handling
pip install pysam

# Real data handling
pip install ont-fast5-api
```

**Hardware Requirements:**
- NVIDIA H100 (current) ✅
- 80GB VRAM ✅
- CUDA 12+ ✅
- Fast SSD storage for genome/reads ✅

**Data Requirements:**
- GRCh38 reference genome ✅ (already have)
- HG002 ONT reads (~50-100 GB)
- minimap2 installed ✅

---

## 7. Files to Create

```
genocache-v4/
├── adaptive_seeding.py         # Adaptive seed extraction
├── seed_chaining.py            # DP-based chaining
├── wfa_gpu_wrapper.py          # WFA GPU integration
├── complete_pipeline.py        # End-to-end pipeline
├── real_data_validation.py     # Real data testing
├── compare_alignments.py       # Compare vs minimap2
├── download_real_data.sh       # Download HG002 data
└── run_complete_benchmark.sh   # Master script
```

---

## 8. Risk Mitigation

**Risk 1: WFA-GPU integration complexity**
- Mitigation: Start with simple test cases
- Fallback: Use CPU WFA initially

**Risk 2: Real data performance differs from synthetic**
- Mitigation: Test on multiple datasets
- Fallback: Retrain model on real data

**Risk 3: Accuracy drops with adaptive seeding**
- Mitigation: Extensive testing and debugging
- Fallback: Keep simple seeding as option

**Risk 4: WFA GPU slower than expected**
- Mitigation: Optimize batching and memory
- Fallback: Use faster alignment algorithm (Edlib)

---

## Next Steps

**IMMEDIATE (after current benchmark):**
1. ✅ Check benchmark results
2. ✅ Validate index quality
3. ✅ Confirm minimap2 comparison

**THEN START:**
1. Implement adaptive seeding
2. Implement seed chaining
3. Integrate WFA GPU
4. Download real ONT data
5. Run complete validation

**ESTIMATED START:** ~4-5 hours from now (after index completes)

---

**Ready to build the complete pipeline! 🚀**

This will take GenoCache from "proof of concept" to "production-ready aligner"!
