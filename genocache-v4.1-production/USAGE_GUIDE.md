# GenoCache v4.1 Usage Guide - minimap2 Compatible

Quick guide for using the minimap2-compatible GenoCache pipeline.

---

## Quick Start

### Basic Usage

```bash
python3 genocache_align_minimap2.py \
    --reads input.fastq \
    --output output.sam \
    --reference GRCh38.fa
```

### With Secondary Alignments

```bash
python3 genocache_align_minimap2.py \
    --reads input.fastq \
    --output output.sam \
    --reference GRCh38.fa \
    --secondary
```

### Custom Model/Index

```bash
python3 genocache_align_minimap2.py \
    --reads input.fastq \
    --output output.sam \
    --reference GRCh38.fa \
    --model models/custom_model.pt \
    --index indexes/custom_index.index \
    --metadata indexes/custom_index.metadata.pkl \
    --secondary
```

---

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--reads` | Input FASTQ file (required) | - |
| `--output` | Output SAM file (required) | - |
| `--reference` | Reference genome FASTA (required) | - |
| `--model` | Model checkpoint path | `models/genocache_model.pt` |
| `--index` | FAISS index path | `indexes/genocache_v4_production.index` |
| `--metadata` | Index metadata path | `indexes/genocache_v4_production.metadata.pkl` |
| `--secondary` | Output secondary alignments | Off |
| `--device` | Device (cpu/cuda) | `cpu` |

---

## Output Format

### SAM Header

```
@HD  VN:1.6  SO:unsorted
@SQ  SN:chr1  LN:248956422
@SQ  SN:chr2  LN:242193529
...
@PG  ID:genocache  PN:genocache  VN:4.1  CL:genocache_align.py
```

### Primary Alignment

```
read_001  0  chr1  12345  60  100M  *  0  0  ACGT...  IIII...  NM:i:2  AS:i:950  ms:i:950  tp:A:P  cm:i:5  s1:i:25  s2:i:800  de:f:0.0132  nn:i:0  rl:i:0
```

### Secondary Alignment

```
read_001  256  chr2  67890  30  98M2I  *  0  0  ACGT...  IIII...  NM:i:5  AS:i:850  ms:i:850  tp:A:S  cm:i:4  s1:i:20  s2:i:0  de:f:0.0245  nn:i:0  rl:i:0
```

### Unmapped Read

```
read_002  4  *  0  0  *  *  0  0  ACGT...  IIII...
```

---

## SAM Tags

All minimap2-compatible tags are included:

| Tag | Type | Description | Example |
|-----|------|-------------|---------|
| `NM` | i | Edit distance | `NM:i:2` |
| `AS` | i | Alignment score | `AS:i:950` |
| `ms` | i | DP alignment score | `ms:i:950` |
| `tp` | A | Type (P/S) | `tp:A:P` |
| `cm` | i | Number of seeds | `cm:i:5` |
| `s1` | i | Best chain score | `s1:i:25` |
| `s2` | i | Second-best chain score | `s2:i:800` |
| `de` | f | Sequence divergence | `de:f:0.0132` |
| `nn` | i | Ambiguous bases | `nn:i:0` |
| `rl` | i | Repetitive seed length | `rl:i:0` |

### MAPQ Values

- **60:** Unique mapping (score ratio > 10)
- **50:** Very confident (score ratio > 5)
- **40:** Confident (score ratio > 3)
- **30:** Good (score ratio > 2)
- **20:** Ambiguous (score ratio > 1.5)
- **10:** Very ambiguous (score ratio > 1.2)
- **0:** Multi-mapping (score ratio ≤ 1.2)

---

## Examples

### Example 1: Basic Alignment

```bash
# Align reads to reference
python3 genocache_align_minimap2.py \
    --reads sample.fastq \
    --output sample.sam \
    --reference GRCh38.fa

# Output:
# Total reads:       100
# Mapped:            73 (73.0%)
# Unmapped:          27 (27.0%)
# Time:              50.9s
# Speed:             1.97 reads/sec
```

### Example 2: With Secondary Alignments

```bash
# Include multi-mapping reads
python3 genocache_align_minimap2.py \
    --reads sample.fastq \
    --output sample.sam \
    --reference GRCh38.fa \
    --secondary

# Output:
# Total reads:       100
# Mapped:            73 (73.0%)
# Unmapped:          27 (27.0%)
# Secondary alns:    31
# Time:              50.9s
# Speed:             1.97 reads/sec
```

### Example 3: Compare with minimap2

```bash
# Run GenoCache
python3 genocache_align_minimap2.py \
    --reads sample.fastq \
    --output genocache.sam \
    --reference GRCh38.fa \
    --secondary

# Run minimap2
minimap2 -ax map-ont GRCh38.fa sample.fastq > minimap2.sam

# Compare
python3 validation/scripts/06_compare_with_minimap2.py
```

---

## Downstream Tool Compatibility

The output SAM format is compatible with:

- ✅ **samtools** - Sort, index, convert to BAM/CRAM
- ✅ **picard** - Mark duplicates, validate SAM
- ✅ **GATK** - Variant calling pipelines
- ✅ **IGV** - Visualization
- ✅ **bedtools** - Coverage analysis
- ✅ **sambamba** - Fast BAM processing
- ✅ **All standard SAM/BAM tools**

### Example: Convert to BAM

```bash
# GenoCache SAM to sorted BAM
samtools view -bS genocache.sam | samtools sort -o genocache.sorted.bam
samtools index genocache.sorted.bam

# Visualize in IGV
igv genocache.sorted.bam
```

---

## Performance Tips

### Speed Optimization

1. **Use GPU if available:**
   ```bash
   --device cuda
   ```

2. **Process in batches:**
   ```bash
   # Split large FASTQ
   split -l 40000 large.fastq batch_
   
   # Process each batch
   for batch in batch_*; do
       python3 genocache_align_minimap2.py --reads $batch --output ${batch}.sam --reference GRCh38.fa
   done
   
   # Merge SAM files
   samtools merge -h header.sam output.bam batch_*.sam
   ```

3. **Adjust FAISS parameters:**
   - Edit `adaptive_seeding.py`
   - Tune `nprobe` (default: 64)
   - Adjust `top_k` (default: 32)

### Memory Management

- Large references: Load chromosomes on-demand
- Many reads: Process in streaming mode
- Low memory: Reduce batch size

---

## Validation Scripts

### Test on Synthetic Data

```bash
# Generate test data
python3 validation/scripts/01_create_test_data.py

# Run baseline test
python3 validation/scripts/03_run_baseline.py

# Run improved test
python3 validation/scripts/04_run_improved.py

# Compare results
python3 validation/scripts/05_compare_results.py
```

### Compare with minimap2

```bash
# Run both aligners
python3 genocache_align_minimap2.py --reads test.fastq --output genocache.sam --reference GRCh38.fa --secondary
minimap2 -ax map-ont GRCh38.fa test.fastq > minimap2.sam

# Compare
python3 validation/scripts/06_compare_with_minimap2.py
```

---

## Troubleshooting

### Low Mapping Rate

**Problem:** Only 73% of reads mapped

**Causes:**
- GenoCache is more conservative than minimap2
- Lower quality reads filtered out
- Stricter alignment thresholds

**Solutions:**
- Adjust `min_score_threshold` in `extend_phase.py`
- Lower `score_ratio_threshold`
- Use real sequencing data for validation

### Slow Processing

**Problem:** ~2 reads/sec (5× slower than minimap2)

**Causes:**
- Neural encoding overhead
- FAISS search time
- Alignment computation

**Solutions:**
- Use GPU (`--device cuda`)
- Optimize FAISS parameters
- Parallelize with batch processing
- Consider minimap2 for high-throughput

### Missing Tags

**Problem:** SAM tags not present

**Causes:**
- Using old `genocache_align.py` instead of new script

**Solutions:**
- Use `genocache_align_minimap2.py`
- Check SAM output with `samtools view -H`

### Wrong Chromosome

**Problem:** Reads mapping to wrong chromosomes

**Causes:**
- Model training data mismatch
- Need more seeds/candidates
- Reference genome issues

**Solutions:**
- Increase `top_k` in seeding
- Lower alignment threshold
- Validate with real data

---

## Known Limitations

1. **Mapping Rate:** 73% vs 100% (minimap2)
   - More conservative thresholds
   - May miss some valid alignments

2. **Speed:** 2 vs 10 reads/sec (minimap2)
   - Neural encoding overhead
   - Not optimized for production scale

3. **Validation:** Synthetic data only
   - Need real GIAB data testing
   - Unknown performance on real reads

4. **Features:** Some minimap2 features missing
   - No supplementary alignments
   - No MD tag
   - No splice alignment

---

## Best Practices

### When to Use GenoCache

✅ **Use GenoCache for:**
- Research applications
- High-accuracy requirements
- Neural alignment experiments
- Small to medium datasets

❌ **Use minimap2 for:**
- Production pipelines
- High-throughput sequencing
- Real-time processing
- Maximum sensitivity needed

### Data Preparation

1. **Quality Control:**
   ```bash
   # Filter low-quality reads first
   seqtk seq -q 20 input.fastq > filtered.fastq
   ```

2. **Reference Genome:**
   - Use same reference as model training
   - Index must match reference version

3. **Read Format:**
   - FASTQ format required
   - Supports gzipped files (.fastq.gz)

---

## Support

### Documentation

- `COMPLETE_VALIDATION_SUMMARY.md` - Full validation results
- `PHASE2_COMPLETE_SUMMARY.md` - minimap2 compatibility
- `CHANGES_LOG.md` - All code changes

### Scripts

- `genocache_align_minimap2.py` - Main pipeline
- `validation/scripts/` - Test and validation scripts

### Issues

- Check alignment logs
- Validate SAM format with samtools
- Compare with minimap2 output

---

## Version History

**v4.1 (2025-11-16):**
- ✅ minimap2-compatible SAM output
- ✅ All required tags (NM, AS, MAPQ, etc.)
- ✅ Secondary alignment support
- ✅ MAPQ calculation (0-60)
- ✅ Enhanced EXTEND phase
- ✅ Comprehensive validation

**v4.0:**
- Initial production version
- Basic SAM output
- WFA2 integration

---

**Quick Reference:**

```bash
# Basic usage
python3 genocache_align_minimap2.py --reads input.fastq --output output.sam --reference GRCh38.fa

# With secondaries
python3 genocache_align_minimap2.py --reads input.fastq --output output.sam --reference GRCh38.fa --secondary

# Convert to BAM
samtools view -bS output.sam | samtools sort -o output.sorted.bam
samtools index output.sorted.bam
```

Ready to align! 🧬
