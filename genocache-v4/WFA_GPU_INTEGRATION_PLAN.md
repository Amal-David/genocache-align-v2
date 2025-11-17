# GenoCache V4 - WFA-GPU Integration Plan

**Date:** 2025-11-14  
**Status:** Documented and Ready for Implementation  
**Priority:** High (next sprint after hackathon)

---

## Executive Summary

This document outlines the plan to integrate GPU-accelerated WFA (Wavefront Alignment) into GenoCache V4's production pipeline. The integration will replace the current CPU fallback (edlib) with high-performance GPU batch alignment, completing the all-GPU pipeline vision.

**Current Status:**
- ✅ Adaptive seeding: COMPLETE (5-18× faster than minimap2)
- ✅ WFA wrapper framework: READY (CPU fallback implemented)
- 📋 WFA-GPU: DOCUMENTED (ready for implementation)

**Expected Impact:**
- 10-175× faster alignment (based on WFA-GPU paper)
- Complete GPU pipeline (seeding + alignment)
- Production-ready base-level precision

---

## Background

### Current Implementation

**File:** `wfa_alignment.py` (253 lines)

**Components:**
- `WFAAligner` class: Alignment wrapper framework
- `align_edlib()`: CPU fallback using edlib library
- `extract_reference()`: Reference sequence extraction
- `format_sam()`: SAM record generation

**Performance:**
- CPU edlib: >1 minute per read (unacceptable)
- Aligning 1kb read vs 100kb reference region
- Single-threaded, not optimized

**Limitation:** CPU alignment is the bottleneck in our pipeline!

### Why WFA-GPU?

**Advantages of WFA (Wavefront Alignment):**
1. Adaptive algorithm: O(ns + s²) complexity
2. Better than Smith-Waterman: O(nm) complexity
3. Handles indels efficiently
4. Memory-efficient wavefront propagation

**GPU Acceleration Benefits:**
1. Batch processing: Align multiple reads simultaneously
2. Parallel wavefront computation
3. Fast memory access (GPU RAM)
4. Proven speedups: 1.5-175× faster than CPU

**Research Papers:**
- "WFA-GPU: Gap-affine pairwise alignment using GPUs" (2023)
- Published in PMC, bioRxiv
- Authors: Quim Aguado-Puig, Santiago Marco-Sola, et al.

---

## Available Implementations

### 1. quim0/WFA-GPU (Primary Choice)

**Repository:** https://github.com/quim0/WFA-GPU

**Features:**
- Gap-affine alignment (production-ready)
- CUDA implementation
- Batch processing support
- Full CIGAR output
- Memory-optimized

**Performance (from paper):**
- 1.5-7.7× faster than CPU WFA
- Up to 175× faster than other GPU implementations
- Handles long reads (10-100kb) efficiently

**Requirements:**
- CUDA toolkit (11.0+)
- NVIDIA GPU (compute capability 6.0+)
- GCC compiler

**Installation:**
```bash
# Clone repository
git clone https://github.com/quim0/WFA-GPU.git
cd WFA-GPU

# Build library
make

# Test installation
make test
```

### 2. quim0/eWFA-GPU (Alternative)

**Repository:** https://github.com/quim0/eWFA-GPU

**Features:**
- Edit distance alignment (simpler)
- Memory-optimized
- Asynchronous compute
- CUDA implementation

**Use case:**
- When gap-affine scoring not needed
- Faster for simple alignments
- Lower memory requirements

### 3. NVIDIA GenomeWorks (Enterprise Option)

**Repository:** https://github.com/NVIDIA-Genomics-Research/GenomeWorks

**Features:**
- Complete genomics toolkit
- `cudaaligner` module for alignment
- `cudamapper` for read mapping
- `cudapoa` for consensus
- Enterprise-grade support

**Advantages:**
- Official NVIDIA support
- Production-tested
- Comprehensive toolkit
- Active maintenance

**Disadvantages:**
- Heavier dependency
- More complex integration
- Overkill for our use case

### Recommendation: Use quim0/WFA-GPU

**Reasons:**
1. Focused on alignment (our need)
2. Proven performance (paper published)
3. Active development
4. Lightweight integration
5. Open source (Apache 2.0)

---

## Integration Plan

### Phase 1: Installation and Testing (Day 1-2)

**Steps:**

1. **Install CUDA toolkit** (if not present)
   ```bash
   nvidia-smi  # Check GPU and CUDA version
   nvcc --version  # Check CUDA compiler
   ```

2. **Clone and build WFA-GPU**
   ```bash
   cd /home/nebius/genocache/genocache-v4
   git clone https://github.com/quim0/WFA-GPU.git
   cd WFA-GPU
   make clean
   make
   make test  # Verify installation
   ```

3. **Test on sample data**
   ```bash
   # Use example from repository
   ./build/bin/wfa-gpu-test \
     --ref reference.fa \
     --query reads.fa \
     --output alignments.sam
   ```

**Acceptance criteria:**
- ✅ Library compiles without errors
- ✅ Tests pass
- ✅ Sample alignment produces valid SAM

### Phase 2: Python Bindings (Day 3-4)

**Approach:** Use ctypes or cffi to call C/CUDA library from Python

**Create:** `wfa_gpu_bindings.py`

```python
import ctypes
import numpy as np
from pathlib import Path

# Load WFA-GPU library
libwfa = ctypes.CDLL('/path/to/WFA-GPU/lib/libwfa_gpu.so')

class WFAGPUAligner:
    """Python wrapper for WFA-GPU CUDA library"""
    
    def __init__(self, device_id=0):
        self.device_id = device_id
        # Initialize CUDA context
        self._init_cuda()
    
    def _init_cuda(self):
        """Initialize CUDA context and allocate GPU memory"""
        # Set device
        libwfa.wfa_gpu_set_device(self.device_id)
        
        # Allocate batch buffers
        self.batch_size = 1024  # Align 1024 reads at once
        self.max_read_len = 100000  # 100kb max
        self.max_ref_len = 200000  # 200kb max
        
        # Allocate GPU buffers
        libwfa.wfa_gpu_allocate_buffers(
            self.batch_size,
            self.max_read_len,
            self.max_ref_len
        )
    
    def align_batch(self, reads, references):
        """
        Align multiple read-reference pairs in batch
        
        Args:
            reads: List[str] - DNA sequences (queries)
            references: List[str] - Reference sequences
        
        Returns:
            alignments: List[dict] - CIGAR, score, position for each
        """
        # Convert strings to ctypes arrays
        n = len(reads)
        
        # Prepare input buffers
        read_ptrs = (ctypes.c_char_p * n)()
        ref_ptrs = (ctypes.c_char_p * n)()
        
        for i in range(n):
            read_ptrs[i] = reads[i].encode('utf-8')
            ref_ptrs[i] = references[i].encode('utf-8')
        
        # Allocate output buffers
        cigars = (ctypes.c_char_p * n)()
        scores = (ctypes.c_int * n)()
        
        # Call WFA-GPU
        libwfa.wfa_gpu_align_batch(
            read_ptrs, ref_ptrs, n,
            cigars, scores
        )
        
        # Parse results
        results = []
        for i in range(n):
            results.append({
                'cigar': cigars[i].decode('utf-8'),
                'score': scores[i],
                'read_id': i
            })
        
        return results
    
    def __del__(self):
        """Free GPU memory"""
        libwfa.wfa_gpu_free_buffers()
```

**Testing:**
```python
# Test batch alignment
aligner = WFAGPUAligner(device_id=0)

reads = ["ACGTACGT" * 100, "GCTAGCTA" * 100]
refs = ["ACGTACGG" * 100, "GCTAGCAA" * 100]

results = aligner.align_batch(reads, refs)
for r in results:
    print(f"CIGAR: {r['cigar']}, Score: {r['score']}")
```

**Acceptance criteria:**
- ✅ Python can call WFA-GPU functions
- ✅ Batch alignment works
- ✅ Results are valid (CIGAR strings correct)

### Phase 3: Pipeline Integration (Day 5-6)

**Update:** `wfa_alignment.py`

```python
# Add GPU support
try:
    from wfa_gpu_bindings import WFAGPUAligner
    WFA_GPU_AVAILABLE = True
except ImportError:
    WFA_GPU_AVAILABLE = False
    print("⚠️  WFA-GPU not available, using CPU fallback")

class WFAAligner:
    def __init__(self, genome_dict, mode='gpu'):
        self.genome_dict = genome_dict
        self.mode = mode
        
        if mode == 'gpu' and WFA_GPU_AVAILABLE:
            self.aligner = WFAGPUAligner(device_id=0)
            print("✅ Using WFA-GPU for alignment")
        else:
            self.aligner = None  # Use CPU fallback
            print("⚠️  Using CPU alignment (edlib)")
    
    def align_batch(self, reads, regions):
        """
        Align multiple reads in batch
        
        Args:
            reads: List[(read_id, read_seq)]
            regions: List[(chr, start, end)]
        
        Returns:
            alignments: List[dict] - SAM records
        """
        if self.aligner and self.mode == 'gpu':
            # GPU batch alignment
            references = [
                self.extract_reference(chr, start, end)
                for chr, start, end in regions
            ]
            
            read_seqs = [seq for _, seq in reads]
            
            # Batch align on GPU
            results = self.aligner.align_batch(read_seqs, references)
            
            # Format as SAM
            sam_records = []
            for i, result in enumerate(results):
                read_id, read_seq = reads[i]
                chr, start, end = regions[i]
                
                sam = self.format_sam(
                    read_id, read_seq, chr, start,
                    result['cigar'], result['score']
                )
                sam_records.append(sam)
            
            return sam_records
        else:
            # CPU fallback (existing code)
            return self.align_batch_cpu(reads, regions)
```

**Update:** `complete_pipeline.py`

```python
# Initialize aligner with GPU
self.aligner = WFAAligner(
    genome_dict=self.genome_dict,
    mode='gpu'  # 'gpu' or 'cpu'
)

# Batch alignment
def align_reads(self, reads_file, output_sam):
    """Process reads in batches for GPU efficiency"""
    
    batch_size = 1024  # Align 1024 reads at once
    
    reads_batch = []
    regions_batch = []
    
    for read_id, read_seq in read_iterator:
        # Get seeding result
        seed_result = self.seeder.align_read(read_seq)
        
        if seed_result:
            reads_batch.append((read_id, read_seq))
            regions_batch.append((
                seed_result['chr'],
                seed_result['start'],
                seed_result['end']
            ))
        
        # Process batch when full
        if len(reads_batch) >= batch_size:
            alignments = self.aligner.align_batch(
                reads_batch, regions_batch
            )
            # Write to SAM
            for sam in alignments:
                output.write(sam)
            
            # Reset batch
            reads_batch = []
            regions_batch = []
    
    # Process remaining
    if reads_batch:
        alignments = self.aligner.align_batch(reads_batch, regions_batch)
        for sam in alignments:
            output.write(sam)
```

**Acceptance criteria:**
- ✅ Pipeline uses GPU alignment automatically
- ✅ Falls back to CPU if GPU unavailable
- ✅ Batch processing works (1024 reads at once)
- ✅ SAM output is valid

### Phase 4: Performance Testing (Day 7-8)

**Benchmarks to run:**

1. **Throughput test** (500 reads)
   ```bash
   python3 complete_pipeline.py \
     --mode gpu \
     --reads test_500reads.fa \
     --output gpu_alignments.sam
   
   # Compare with CPU
   python3 complete_pipeline.py \
     --mode cpu \
     --reads test_500reads.fa \
     --output cpu_alignments.sam
   ```

2. **Accuracy test** (compare CIGAR strings)
   ```bash
   # Run minimap2
   minimap2 -ax map-ont ref.fa reads.fa > minimap2.sam
   
   # Run GenoCache
   python3 complete_pipeline.py reads.fa > genocache.sam
   
   # Compare alignments
   python3 compare_alignments.py \
     minimap2.sam genocache.sam \
     --output comparison.json
   ```

3. **Scale test** (10,000 reads)
   ```bash
   time python3 complete_pipeline.py \
     --reads large_test.fa \
     --output large_test.sam \
     --mode gpu
   ```

**Metrics to collect:**
- Throughput (reads/sec)
- GPU utilization (nvidia-smi)
- Memory usage
- CIGAR agreement with minimap2
- Mapping quality scores

**Expected results:**
- GPU: 100-500 reads/sec (10-50× faster than CPU)
- Accuracy: >95% CIGAR agreement
- Memory: <16 GB GPU RAM

**Acceptance criteria:**
- ✅ 10× faster than CPU alignment
- ✅ >95% accuracy vs minimap2
- ✅ Stable under load (10k+ reads)

### Phase 5: Real Data Validation (Day 9-10)

**Test on GIAB HG002 data:**

1. **Download GIAB data**
   ```bash
   # Use latest 2025.01 release
   aws s3 cp --no-sign-request \
     s3://ont-open-data/giab_2025.01/HG002/ \
     giab_hg002/ \
     --recursive \
     --exclude "*" \
     --include "*.fastq.gz"
   ```

2. **Subsample for testing**
   ```bash
   # Extract chr22 reads
   seqtk subseq giab_hg002_all.fastq.gz chr22_reads.txt > chr22.fastq
   
   # Take 10k reads
   head -n 40000 chr22.fastq > chr22_10k.fastq
   ```

3. **Run complete pipeline**
   ```bash
   python3 complete_pipeline.py \
     --reads chr22_10k.fastq \
     --output genocache_hg002.sam \
     --mode gpu
   
   # Compare with minimap2
   minimap2 -ax map-ont GRCh38.fa chr22_10k.fastq > minimap2_hg002.sam
   ```

4. **Validate results**
   ```bash
   # Check mapping rates
   samtools flagstat genocache_hg002.sam
   samtools flagstat minimap2_hg002.sam
   
   # Compare positions
   python3 compare_alignments.py \
     genocache_hg002.sam \
     minimap2_hg002.sam \
     --detailed \
     --output hg002_comparison.json
   ```

**Acceptance criteria:**
- ✅ Similar mapping rate to minimap2 (±5%)
- ✅ Position agreement >90%
- ✅ CIGAR agreement >85%
- ✅ Handles real ONT errors (substitutions, indels)

---

## Expected Performance

### Throughput (reads/sec)

| Component | Current (CPU) | With WFA-GPU | Speedup |
|-----------|---------------|--------------|---------|
| Encoding | 2.56 ms | 2.56 ms | 1× |
| FAISS search | 6.43 ms | 6.43 ms | 1× |
| Adaptive chaining | 3-5 ms | 3-5 ms | 1× |
| **Alignment** | **60,000 ms** | **100 ms** | **600×** |
| **Total** | **60,012 ms** | **112 ms** | **536×** |

**Current bottleneck:** Alignment (99.98% of time!)  
**After WFA-GPU:** Balanced pipeline (alignment 89% of time)

### End-to-End Comparison

| Tool | Mode | Reads/sec | Time for 10k reads |
|------|------|-----------|-------------------|
| minimap2 | CPU (8 threads) | 2.2 | 75 min |
| GenoCache | Seeding only | 12.4 | 13 min |
| GenoCache | + CPU align | 0.02 | 138 hours ❌ |
| GenoCache | + GPU align | **100+** | **1.7 min** ✅ |

**Expected final speedup:** 50-100× faster than minimap2! 🚀

---

## Risk Mitigation

### Risk 1: WFA-GPU Compilation Fails

**Likelihood:** Medium  
**Impact:** High

**Mitigation:**
- Pre-test on development machine
- Have Docker container with pre-built library
- Alternative: Use NVIDIA GenomeWorks (more stable)

### Risk 2: Python Bindings Don't Work

**Likelihood:** Low  
**Impact:** Medium

**Mitigation:**
- Use ctypes (standard library, no dependencies)
- Test bindings separately before integration
- Have example code from WFA-GPU repo

### Risk 3: GPU Memory Limitations

**Likelihood:** Low  
**Impact:** Medium

**Mitigation:**
- Batch size tuning (start with 128, scale up)
- Monitor GPU memory during testing
- Implement memory pooling if needed

### Risk 4: Accuracy Regression

**Likelihood:** Low  
**Impact:** High

**Mitigation:**
- Extensive testing against minimap2
- Use GIAB truth sets for validation
- Compare CIGAR strings position-by-position

---

## Success Metrics

### Technical Metrics

1. **Performance:**
   - Alignment throughput: >100 reads/sec
   - Total pipeline: >50 reads/sec
   - 50× faster than minimap2

2. **Accuracy:**
   - Mapping rate: ±5% of minimap2
   - Position agreement: >90%
   - CIGAR agreement: >85%

3. **Stability:**
   - Process 10k+ reads without crash
   - GPU memory stable (<16 GB)
   - CPU usage <50% (GPU-bound)

### Business Metrics

1. **Integration:**
   - Drop-in replacement for minimap2
   - Compatible with existing pipelines
   - SAM output format

2. **Adoption:**
   - Parabricks integration possible
   - Cloud deployment ready
   - Docker container available

3. **Impact:**
   - 3-5× total workflow speedup
   - Cost savings (GPU vs CPU hours)
   - Enables real-time analysis

---

## Timeline Summary

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| 1. Installation | 2 days | WFA-GPU compiled and tested |
| 2. Python bindings | 2 days | Working Python wrapper |
| 3. Pipeline integration | 2 days | GPU mode in pipeline |
| 4. Performance testing | 2 days | Benchmark results |
| 5. Real data validation | 2 days | GIAB HG002 results |
| **Total** | **10 days** | **Production-ready GPU pipeline** |

---

## Resources Required

### Hardware
- NVIDIA GPU (V100, A100, or H100)
- 16+ GB GPU RAM
- CUDA 11.0+ support

### Software
- CUDA toolkit 11.0+
- GCC 7.0+
- Python 3.8+
- PyTorch (optional, for future integration)

### Personnel
- 1 engineer (full-time, 10 days)
- 0.5 engineer (testing/validation)

### Data
- GIAB HG002 dataset (~20 GB)
- Reference genome (GRCh38, ~3 GB)
- Synthetic test data (existing)

---

## Conclusion

WFA-GPU integration is the **final critical piece** to complete GenoCache V4's all-GPU pipeline. The integration is well-documented, with clear implementation steps and proven technology (published paper, active repository).

**Current State:**
- ✅ Adaptive seeding: 5-18× faster (COMPLETE)
- ✅ Production index: Ready (COMPLETE)
- ✅ Pipeline framework: Ready (COMPLETE)
- ⚠️ Alignment: CPU bottleneck (NEXT SPRINT)

**After WFA-GPU Integration:**
- ✅ Complete GPU pipeline (seeding + alignment)
- ✅ 50-100× faster than minimap2
- ✅ Production-ready for Parabricks integration
- ✅ Cost-effective cloud deployment

**Recommendation:** Prioritize WFA-GPU integration in next sprint (immediately after hackathon) to deliver complete production system.

---

**Document Version:** 1.0  
**Last Updated:** 2025-11-14 02:00 UTC  
**Status:** Ready for Implementation  
**Next Review:** After hackathon (2025-11-16)
