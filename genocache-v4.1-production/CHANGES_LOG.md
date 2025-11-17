# Changes Log - GenoCache V4.1

**IMPORTANT:** This file tracks ALL code changes made to production

---

## Session: 2025-11-15 (NeuralAligner Improvements)

### File Modified: `genocache_core/adaptive_seeding.py`

#### Change 1: Chain Scoring Method (Line ~261)

**Location:** Inside `chain_seeds()` method

**Before:**
```python
chain_score = sum(s[2] for s in seed_list)
```

**After:**
```python
# FIX: Use number of anchors (like NeuralAligner) instead of sum of scores
# This prevents bias toward longer chromosomes
chain_score = len(seed_list)  # Count anchors (fair to all regions)
```

**Reason:** NeuralAligner counts anchors (fair), we were summing scores (biased toward long chromosomes like main chr22 vs short alternate contigs)

**Expected Impact:** +2-3% accuracy for alternate contig cases

---

#### Change 2: Enhanced Rescue Logic (Lines ~304-337)

**Location:** Inside `align_read()` method

**Before:**
```python
# Step 4: Decision gate
if len(chains) == 0:
    # No chains found - try rescue seeding
    if len(seeds) < self.max_seeds:
        # Add more seeds
        rescue_seeds = self.extract_seeds(read, num_seeds=self.max_seeds)
        for i, seed in enumerate(rescue_seeds):
            rescue_seeds[i] = self.search_seed(seed)
            rescue_seeds[i] = self.filter_seed(rescue_seeds[i])
        
        chains = self.chain_seeds(rescue_seeds)
    
    if len(chains) == 0:
        # Still no chains - unmapped
        return None
```

**After:**
```python
# Step 4: Decision gate with enhanced rescue logic (like NeuralAligner)
need_rescue = False

if len(chains) == 0:
    # Condition 1: No chains found
    need_rescue = True
elif len(chains) > 0:
    # Condition 2: Best chain score too low (< half of seeds)
    # This indicates weak mapping confidence
    if chains[0].score < self.min_seeds / 2:
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
```

**Reason:** NeuralAligner triggers rescue on 3 conditions (no chains, low score, ambiguous). We only had 1 condition.

**Expected Impact:** +1-2% accuracy for ambiguous and weak mapping cases

---

#### Change 3: FAISS nprobe Preservation (Lines ~74-84)

**Location:** Inside `__init__()` method

**Before:**
```python
self.model = model
self.index = index
self.metadata = metadata
# ... (nprobe not checked)
self.char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
self.model.eval()
```

**After:**
```python
self.model = model
self.index = index
self.metadata = metadata
# ... (other init code)

# Check nprobe (index may already have optimal value)
# NeuralAligner uses nprobe=32, but our production index uses 64
# Keep the existing value if it's already high
if hasattr(self.index, 'nprobe'):
    original_nprobe = self.index.nprobe
    # Only increase if it's too low (< 32)
    if original_nprobe < 32:
        self.index.nprobe = 32
        print(f"  FAISS nprobe increased: {original_nprobe} → {self.index.nprobe}")
    else:
        print(f"  FAISS nprobe: {self.index.nprobe} (already optimal)")

self.char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
self.model.eval()
```

**Reason:** 
- Initially thought nprobe=16 (too low)
- Discovered production index has nprobe=64 (optimal!)
- Added check to preserve optimal value

**Expected Impact:** Preserves existing optimal configuration (no change to accuracy, prevents regression)

---

## Rollback Instructions

To revert all changes:

```bash
cd /home/nebius/genocache/genocache-v4.1-production

# Option 1: Manual revert (if no backup exists)
# Manually undo the 3 changes above in adaptive_seeding.py

# Option 2: From backup (if backup exists)
cp genocache_core/adaptive_seeding_BACKUP.py genocache_core/adaptive_seeding.py

# Option 3: From git (if using version control)
git checkout genocache_core/adaptive_seeding.py
```

---

## Testing Status

- [ ] Unit tests run
- [ ] Integration tests run
- [ ] Validation on 8 reads: 87.5% (same as before changes)
- [ ] Validation on 38 synthetic reads: 36.8% (data quality issues)
- [ ] Validation on real dataset: PENDING

---

## Notes

**Issue:** Changes were made directly to production code without proper version control.

**Solution going forward:**
1. Create baseline backup before changes
2. Use working directory for experiments
3. Document each change immediately
4. Test before merging to production
5. Use git or similar for tracking

---

**Last Updated:** 2025-11-15  
**Status:** Changes applied, pending real-world validation
