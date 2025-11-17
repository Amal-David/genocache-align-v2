# GenoCache V4.1 - Quick Start Guide

Get up and running with GenoCache in 5 minutes!

---

## Step 1: Install Dependencies

```bash
# Create virtual environment (recommended)
python3 -m venv genocache-env
source genocache-env/bin/activate  # On Windows: genocache-env\Scripts\activate

# Install required packages
pip install torch faiss-cpu biopython parasail-python numpy
```

### Dependency Details

- **torch** - PyTorch for neural network
- **faiss-cpu** - FAISS for vector search (use faiss-gpu if you have CUDA)
- **biopython** - For FASTA/FASTQ parsing
- **parasail-python** - For alignment (SSE/AVX optimized)
- **numpy** - For numerical operations

---

## Step 2: Download Reference Genome

```bash
# Download GRCh38 (if you don't have it)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Uncompress
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa
```

**Note:** Place this in a known location, we'll reference it in the pipeline.

---

## Step 3: Prepare Your Reads

GenoCache accepts FASTA or FASTQ format:

```bash
# Example: Your reads file
your_reads.fastq  # or your_reads.fasta
```

**Supported formats:**
- FASTQ (.fastq, .fq)
- FASTA (.fasta, .fa)
- Compressed versions (.gz) - decompress first

---

## Step 4: Run GenoCache

```bash
# Basic usage
python genocache_align.py \
    --reads your_reads.fastq \
    --output aligned.sam \
    --reference /path/to/GRCh38.fa
```

### Expected Output

```
================================================================================
GenoCache V4.1 Production Pipeline
================================================================================

Loading model from models/genocache_model.pt...
✅ Model loaded
Loading FAISS index from indexes/genocache_v4_production.index...
✅ Index loaded: 91,792,546 vectors
Loading metadata from indexes/genocache_v4_production.metadata.pkl...
✅ Metadata loaded
Loading reference genome from /path/to/GRCh38.fa...
  NC_000001.11: 248,956,422 bp
  NC_000002.12: 242,193,529 bp
  ...
✅ Reference loaded: 379 sequences

Initializing pipeline...
✅ Adaptive seeder initialized
✅ Fast aligner initialized
✅ EXTEND phase initialized

================================================================================
Processing reads
================================================================================

Loaded 1000 reads

  Processed 100/1000 reads...
  Processed 200/1000 reads...
  ...

================================================================================
Alignment Complete
================================================================================
Total reads: 1000
Mapped: 875 (87.5%)
Unmapped: 125 (12.5%)
Time: 500.0s (2.0 reads/sec)
Output: aligned.sam
```

---

## Step 5: Check Results

```bash
# View SAM file
head -n 50 aligned.sam

# Count mapped vs unmapped
grep -v "^@" aligned.sam | awk '{if ($2 == 0 || $2 == 16) print "mapped"; else print "unmapped"}' | sort | uniq -c
```

### SAM Format

GenoCache produces standard SAM format:

```
@HD	VN:1.0	SO:unsorted
@SQ	SN:NC_000001.11	LN:248956422
@SQ	SN:NC_000002.12	LN:242193529
...
@PG	ID:genocache	PN:GenoCache	VN:4.1.0
read_1	0	NC_000022.11	34371979	60	1000M	*	0	0	ATACG...	*
read_2	4	*	0	0	*	*	0	0	CGTAT...	*
...
```

**Fields:**
- Column 1: Read ID
- Column 2: FLAG (0=mapped forward, 4=unmapped, 16=mapped reverse)
- Column 3: Reference chromosome
- Column 4: Position
- Column 5: MAPQ (mapping quality)
- Column 6: CIGAR string
- Column 10: Read sequence

---

## Advanced Usage

### Custom Parameters

```bash
# Adjust EXTEND phase candidates
python genocache_align.py \
    --reads your_reads.fastq \
    --output aligned.sam \
    --reference /path/to/GRCh38.fa \
    --top-k 10  # Test more candidates (default: 5)

# Use custom model/index
python genocache_align.py \
    --reads your_reads.fastq \
    --output aligned.sam \
    --reference /path/to/GRCh38.fa \
    --model /path/to/custom_model.pt \
    --index /path/to/custom.index \
    --metadata /path/to/custom.metadata.pkl
```

### Disable EXTEND (Not Recommended)

```bash
# Use old method (for comparison)
python genocache_align.py \
    --reads your_reads.fastq \
    --output aligned_old.sam \
    --reference /path/to/GRCh38.fa \
    --no-extend
```

**Note:** This will give 37.5% accuracy instead of 87.5%. Only use for debugging.

---

## Validation

### Run Built-in Validation

```bash
# Validate EXTEND fix on test data
python scripts/validate_extend_fix.py
```

This will show:
- OLD method accuracy (37.5%)
- NEW method accuracy (87.5%)
- Alignment scores for each candidate
- Which reads were fixed by EXTEND

---

## Troubleshooting

### Issue: "Model not found"

```
Error: Model file not found: models/genocache_model.pt
```

**Solution:** Check that you're running from the `genocache-v4.1-production/` directory:

```bash
cd /path/to/genocache-v4.1-production
python genocache_align.py --reads ...
```

### Issue: "FAISS index not found"

```
Error: FAISS index not found: indexes/genocache_v4_production.index
```

**Solution:** Make sure you have the complete package with all files. Check:

```bash
ls -lh indexes/
# Should show:
# genocache_v4_production.index (2.1GB)
# genocache_v4_production.metadata.pkl (4.8GB)
```

### Issue: "Reference file not found"

```
Error: Reference file not found: /path/to/GRCh38.fa
```

**Solution:** Download reference genome (see Step 2) and use full path:

```bash
python genocache_align.py \
    --reference /absolute/path/to/GRCh38.fa \
    ...
```

### Issue: "Out of memory"

```
RuntimeError: CUDA out of memory
```

**Solution:** GenoCache runs on CPU by default. If you modified it to use GPU:

```python
# In genocache_align.py, keep model on CPU:
model = model.cpu()  # Not .cuda()
```

### Issue: "Slow performance"

**Expected:** ~1-2 reads/sec with Parasail

**If slower:**
- Check CPU load (other processes?)
- Check disk I/O (reference genome on slow drive?)
- Try smaller batch first to isolate issue

**For faster performance:** Wait for WFA-GPU integration (250× speedup coming soon!)

---

## Next Steps

### Compare with minimap2

```bash
# Run minimap2 on same reads
minimap2 -ax sr /path/to/GRCh38.fa your_reads.fastq > minimap2.sam

# Compare chromosome accuracy
python scripts/compare_sam_files.py \
    --genocache aligned.sam \
    --minimap2 minimap2.sam
```

### Large-scale Testing

```bash
# Test on full dataset
python genocache_align.py \
    --reads large_dataset.fastq \
    --output large_aligned.sam \
    --reference /path/to/GRCh38.fa

# Monitor progress
tail -f genocache.log  # if logging enabled
```

---

## Performance Tips

### 1. Use SSD for Reference

Place GRCh38.fa on SSD for faster access:
- HDD: ~100 MB/s
- SSD: ~500 MB/s
- NVMe: ~3000 MB/s

### 2. Preload Reference

For multiple runs, keep reference in memory (use RAMdisk):

```bash
# Linux
sudo mkdir /mnt/ramdisk
sudo mount -t tmpfs -o size=4G tmpfs /mnt/ramdisk
cp GRCh38.fa /mnt/ramdisk/

# Use ramdisk reference
python genocache_align.py --reference /mnt/ramdisk/GRCh38.fa ...
```

### 3. Batch Processing

For large datasets, split into batches:

```bash
# Split reads
split -l 40000 large_reads.fastq batch_  # 10k reads per batch

# Process each
for batch in batch_*; do
    python genocache_align.py \
        --reads $batch \
        --output ${batch}.sam \
        --reference /path/to/GRCh38.fa
done

# Merge SAM files
# (Keep only one header, concatenate alignments)
```

---

## Getting Help

### Documentation

- `README.md` - Overview and features
- `docs/EXTEND_PHASE_VALIDATION.md` - Technical validation
- `docs/QUICK_START.md` - This guide

### Support

For issues:
1. Check troubleshooting section above
2. Review error messages carefully
3. Check GitHub Issues
4. Contact support email

### Feedback

We want to hear from you!
- Report bugs
- Request features
- Share results
- Contribute code

---

## Summary

You're now ready to use GenoCache V4.1!

**Quick recap:**
1. ✅ Install dependencies (`pip install torch faiss-cpu biopython parasail-python`)
2. ✅ Download reference genome (GRCh38.fa)
3. ✅ Run pipeline (`python genocache_align.py --reads ... --output ... --reference ...`)
4. ✅ Check results (SAM file with CIGAR strings)

**Key features:**
- 87.5% chromosome accuracy (vs 37.5% without EXTEND)
- Standard SAM output
- Production-ready
- WFA-GPU acceleration coming soon!

---

**Happy Aligning!** 🧬🚀
