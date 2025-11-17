# Crash Prevention & Recovery Guide

**All improvements implemented in `encode_grch38_robust.py`**

---

## 🛡️ Anti-Crash Features Implemented

### 1. **Per-Chromosome Checkpointing** ✅
```python
# Each chromosome saved immediately after encoding
encoding_checkpoints/
  chr_0_NC_000001.11.npz  ← Chr 1 complete
  chr_1_NC_000002.12.npz  ← Chr 2 complete
  ...
```

**Benefits:**
- No work lost if crash occurs
- Resume from last completed chromosome
- Can inspect progress anytime

### 2. **Comprehensive Error Logging** ✅
```python
# Tracks two categories:

skipped_chromosomes = []  # Expected skips (scaffolds, etc.)
  - "NT_187361.1 (not main chromosome)"
  - "scaffold_123 (too short: 1000 bp)"
  - "chrUn_KI270442v1 (too many Ns)"

failed_chromosomes = []   # Unexpected failures
  - "NC_000001.11 (CUDA out of memory)"
  - "NC_000015.10 (checkpoint save failed)"
  - "NC_000023.11 (encoding error: ...)"
```

**Output:**
- Printed in console at end
- Saved to `encoding_failures.log`
- Helps debug issues

### 3. **Memory Management** ✅
```python
# After each chromosome:
del chr_vectors, chr_positions, chr_chromosomes, windows
gc.collect()  # Python garbage collection
torch.cuda.empty_cache()  # GPU memory cleanup

# During encoding:
if i % (BATCH * 10) == 0:
    torch.cuda.empty_cache()  # Clear every 10 batches
```

### 4. **Safe File Operations** ✅
```python
# Checkpoint saves use compression
np.savez_compressed(checkpoint_file, ...)

# Final saves have try-except
try:
    np.save(OUT_VECTORS, vectors)
except Exception as e:
    print(f"Save failed: {e}")
    # Files remain in checkpoints/
```

### 5. **Resume Capability** ✅
```python
# Automatic checkpoint detection
if checkpoints exist:
    load_checkpoints()
    resume from last chromosome
else:
    start fresh
```

---

## 🔧 Recovery Procedures

### If Encoding Crashes Mid-Run

**Option A: Resume Automatically**
```bash
# Just re-run the same script
.venv/bin/python encode_grch38_robust.py

# It will:
# - Detect existing checkpoints
# - Skip completed chromosomes
# - Continue where it left off
```

**Option B: Assemble from Checkpoints**
```bash
# If encoding is complete but final save failed
.venv/bin/python resume_from_checkpoints.py

# This will:
# - Load all checkpoint files
# - Assemble into final .npy files
# - No re-encoding needed
```

**Option C: Check What's Done**
```bash
# Count completed chromosomes
ls encoding_checkpoints/*.npz | wc -l

# List them
ls -lh encoding_checkpoints/

# Check total vectors so far
python3 << 'EOF'
import glob, numpy as np
files = sorted(glob.glob('encoding_checkpoints/*.npz'))
total = sum(len(np.load(f)['vectors']) for f in files)
print(f"Total vectors in checkpoints: {total:,}")
EOF
```

---

## 📋 Monitoring Commands

### Real-Time Progress
```bash
# Watch checkpoint count increase
watch -n 10 'echo "Checkpoints: $(ls encoding_checkpoints/*.npz 2>/dev/null | wc -l)/25"; tail -5 grch38_robust.log'

# Follow log
tail -f grch38_robust.log

# Check process
ps aux | grep encode_grch38_robust
```

### Check Health
```bash
# GPU memory
nvidia-smi

# Disk space
df -h /home/nebius/genocache

# RAM usage
free -h

# Process memory
ps -p $(cat grch38_robust.pid) -o pid,vsz,rss,cmd
```

---

## ⚠️ Common Issues & Solutions

### Issue 1: Out of Disk Space
```bash
# Symptoms: "No space left on device"

# Check space
df -h

# Quick fix: Clean old files
rm -f *.log.old backup_*/ train.log
rm -f faiss_ivfflat_*.idx faiss_ivfpq_M4_*.idx

# Extreme: Remove old vectors
rm -f ref_vectors.npy ref_positions.npy
# (Keep ref_vectors_improved.npy!)
```

### Issue 2: Out of Memory (RAM)
```bash
# Symptoms: "Killed" or "MemoryError"

# Check memory
free -h

# Solutions:
1. Reduce BATCH size in script:
   Edit encode_grch38_robust.py
   Change: BATCH = 2048 → BATCH = 1024

2. Enable swap (if available)

3. Process fewer chromosomes at once
```

### Issue 3: GPU Out of Memory
```bash
# Symptoms: "CUDA out of memory"

# Solutions:
1. Reduce BATCH size:
   BATCH = 2048 → BATCH = 512

2. Clear cache more frequently:
   # Already done in robust script

3. Use CPU if desperate:
   DEVICE = "cpu"  # Slower but works
```

### Issue 4: Checkpoint Save Fails
```bash
# Symptoms: "checkpoint save failed" in log

# Check:
ls -lh encoding_checkpoints/
df -h  # Disk space

# The script will:
- Log the failure
- Skip that chromosome
- Continue with others
- You can retry failed ones later
```

---

## 🎯 Current Robust Encoding Status

**Running**: PID 183840  
**Script**: `encode_grch38_robust.py`  
**Log**: `grch38_robust.log`

**Improvements over original:**
1. ✅ Checkpoints every chromosome (can resume)
2. ✅ Tracks all skipped chromosomes (with reasons)
3. ✅ Tracks all failed chromosomes (with errors)
4. ✅ Saves failure log automatically
5. ✅ Better memory management
6. ✅ GPU cache clearing
7. ✅ Can resume if interrupted
8. ✅ Assembles from checkpoints if final save fails

**Maximum possible data loss:** One chromosome (the one currently encoding)  
**Recovery time:** 0 minutes (just resume)

---

## 📊 What Gets Logged

### Console Output (and grch38_robust.log):
```
Processing NC_000001.11: 248,956,422 bp
  7,202,183 windows extracted
  Encoding... done (7,202,183 vectors)
  ✓ Checkpoint saved: encoding_checkpoints/chr_0_NC_000001.11.npz
  Progress: 0.25 Gbp @ 0.7 Mbp/s
```

### End Summary:
```
Total chromosomes: 25
Total vectors: 89,123,456

Skipped chromosomes (143):
  - NT_113878.1 (not main chromosome)
  - NW_003315944.1 (not main chromosome)
  ...

⚠️  Failed chromosomes (0):
  (none - all succeeded!)
```

### If Failures Occur (encoding_failures.log):
```
Failed Chromosomes:
NC_000015.10 (encoding error: CUDA out of memory at batch 500)
NC_000023.11 (checkpoint save failed: disk full)
```

---

## ✅ Verification Checklist

### After Encoding Completes

**Check 1: Count checkpoints**
```bash
ls encoding_checkpoints/*.npz | wc -l
# Should be: 25 (main chromosomes)
```

**Check 2: Check log for failures**
```bash
grep "Failed chromosomes" grch38_robust.log
# Should show: (0) if all succeeded
```

**Check 3: Verify final files exist**
```bash
ls -lh grch38_*.npy
# Should show 3 files: vectors, positions, chromosomes
```

**Check 4: Verify file sizes**
```bash
# Vectors should be 30-50 GB
# Positions should be 30-100 MB  
# Chromosomes should be 30-100 MB
```

**Check 5: Quick sanity check**
```python
import numpy as np
v = np.load('grch38_vectors.npy', mmap_mode='r')
p = np.load('grch38_positions.npy')
c = np.load('grch38_chromosomes.npy')

print(f"Vectors: {v.shape}")
print(f"Positions: {len(p):,}")
print(f"Chromosomes: {len(np.unique(c))}")
print(f"All lengths match: {len(v) == len(p) == len(c)}")
```

---

## 🚀 Quick Start Recovery

### If Process Dies
```bash
# Check what's saved
ls encoding_checkpoints/*.npz | wc -l

# If less than 25: Resume encoding
.venv/bin/python encode_grch38_robust.py

# If 25 checkpoints exist: Assemble finals
.venv/bin/python resume_from_checkpoints.py
```

### Force Clean Restart
```bash
# Remove all checkpoints
rm -rf encoding_checkpoints/

# Remove logs
rm grch38_robust.log grch38_robust.pid

# Start fresh
.venv/bin/python encode_grch38_robust.py
```

---

**Current Status:** Robust encoding running with full crash protection! 🛡️  
**All failures will be logged with reasons** ✓  
**No work will be lost** ✓  
**Can resume anytime** ✓

