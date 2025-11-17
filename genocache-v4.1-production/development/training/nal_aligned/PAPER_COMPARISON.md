# NAL Paper vs Our Implementation - Technical Comparison

## CRITICAL DIFFERENCES ⚠️

### 1. **Seed Numbers (Section A.5)**

**Paper:**
```
"In our implementation, we start with 7 seeds per read and 
increase to 13 seeds in the second round."
```

**Our Implementation:**
```python
# align_nal.py, line 105-106
anchors_iter1 = self.seeder.get_anchors(
    read_seq, num_seeds=5, K=self.K  # ❌ WRONG: Should be 7
)
...
anchors_iter2 = self.seeder.get_anchors(
    read_seq, num_seeds=16, K=self.K  # ❌ WRONG: Should be 13
)
```

**Status:** ❌ **MISMATCH** - We use 5→16, paper uses **7→13**

---

### 2. **Chaining Implementation (Equation 2)**

**Paper Equation 2:**
```
y_chain^+ = argmax_y [ Σ_{i=1}^n Σ_{j=1}^m I(|y_ij - x_i - y| ≤ C) ]
y_chain^- = argmax_y [ Σ_{i=1}^n Σ_{j=1}^m I(|y_ij + x_i - L_read + L_seed - y| ≤ C) ]
```

**Our Implementation:**
```python
# chaining_nal.py, lines 143-150
if strand == '+':
    # Positive strand: y = y_ij - x_i
    expected_y = y_ij - x_i  # ✓ CORRECT
else:
    # Negative strand: y = y_ij + x_i - read_len + seed_len
    expected_y = y_ij + x_i - read_len + seed_len  # ✓ CORRECT
```

**Status:** ✅ **CORRECT** - Matches Equation 2 exactly

---

### 3. **Chaining Score (Section 3.4)**

**Paper:**
```
"The optimal chain is defined as the one that contains the 
largest number of anchors rather than the longest span."
```

**Our Implementation:**
```python
# chaining_nal.py, line 168
score = len(chain_anchors)  # ✓ CORRECT - count anchors
```

**Status:** ✅ **CORRECT**

---

### 4. **Indexing Parameters (Section A.3)**

| Parameter | Paper Spec | Our Implementation | Status |
|-----------|------------|-------------------|--------|
| Index type | IVFPQ | IVFPQ | ✅ |
| nlist | sqrt(N) | sqrt(81.6M) = 9035 | ✅ |
| Stride | ≤ seed_len/8 (≤64 for 512) | 32 | ✅ |
| nprobe | 8-32 | 8 | ✅ |
| K neighbors | 32 | 32 | ✅ |
| Distance | Inner product | Inner product | ✅ |
| PQ compression | PQ16×8 | PQ16×8 | ✅ |

**Status:** ✅ **ALL CORRECT**

---

### 5. **Seed Length**

**Paper (Section 3.3):**
```
"typically 256 or 512"
```

**Paper (Section A.2):**
```
"L ranges from 64 to 1024"
```

**Paper Experiments (Table 4):**
```
Tests with L_seed = 128, 256, 384, 512
```

**Our Implementation:**
```python
seed_len = 512  # ✓ CORRECT - within paper's range
```

**Status:** ✅ **CORRECT**

---

### 6. **Tolerance C (Section 3.4)**

**Paper:**
```
"C is the allowed positional tolerance"
(Algorithm 1 uses C_b for bias tolerance)
```

**Typical value:** Not explicitly stated in main text, but implied ~1000bp from context

**Our Implementation:**
```python
tolerance = 1000  # ✓ REASONABLE - typical for long reads
```

**Status:** ⚠️ **REASONABLE** - Not explicitly specified in paper

---

### 7. **Top-K Chains**

**Paper (Section 3.4):**
```
"In practice, we retain the top-K chains with the most anchors 
for the next stage of alignment."
```

**Our Implementation:**
```python
top_k = 5  # Default in chaining_nal.py
```

**Status:** ⚠️ **NOT SPECIFIED** - Paper doesn't give exact K value

---

### 8. **WFA Alignment (Section 3.6)**

**Paper:**
```
"we extract from the reference a segment of length 1.002 × L_read"
```

**Our Implementation:**
```python
# align_nal.py, line 268
ref_len = int(read_len * 1.002)  # ✅ CORRECT
```

**Status:** ✅ **CORRECT** (but WFA is placeholder)

---

### 9. **Training Parameters (Section A.2)**

| Parameter | Paper | Our Implementation | Status |
|-----------|-------|-------------------|--------|
| Dimension | D=128 | 128 | ✅ |
| Batch size | 8192 | 8192 | ✅ |
| Temperature τ | Not specified | 0.07 | ⚠️ |
| Optimizer | AdamW | AdamW | ✅ |
| Learning rate | 1e-3 | 1e-3 | ✅ |
| Weight decay | 0.1 | 0.1 | ✅ |
| Error rate | U[0.01, 0.1] | U[0.01, 0.1] | ✅ |
| Shift | ±L_seed/10 | ±51bp (512/10) | ✅ |

**Status:** ✅ **MOSTLY CORRECT**

---

## SUMMARY

### ❌ Critical Issue: Seed Numbers

**Paper uses 7→13, we use 5→16**

This is from Section A.5:
```
"For speed experiments, we use the fixed-length and identity dataset. 
Indexes are identical to those in the accuracy experiments, except that 
the number of seeds per iteration is set to 5 and 11."

"For alignment accuracy experiments... We use 7 seeds in the first 
iteration and 13 seeds in the second."
```

**Impact:** Different seed numbers may affect:
- Mapping accuracy (paper achieves 99.6-100%)
- Speed (more seeds = slower)
- Rescue frequency

### ✅ Core Algorithm: Correct

- Chaining equation matches exactly
- Score = anchor count (not span)
- Index parameters all correct
- Training protocol correct

### ⚠️ Minor Differences

- Tolerance C not explicitly specified (we use 1000bp)
- Top-K chains not specified (we use K=5)
- Temperature τ not specified (we use 0.07)

---

## RECOMMENDATION

### **Fix seed numbers to match paper:**

```python
# Change in align_nal.py:

# ITERATION 1: 7 seeds (not 5)
anchors_iter1 = self.seeder.get_anchors(
    read_seq, num_seeds=7, K=self.K  # ✓ Match paper
)

# ITERATION 2: 13 seeds (not 16)
if needs_rescue:
    anchors_iter2 = self.seeder.get_anchors(
        read_seq, num_seeds=13, K=self.K  # ✓ Match paper
    )
```

### **Why this matters:**

1. **Accuracy comparison:** Paper reports 99.6-100% with 7→13 seeds
2. **Reproducibility:** Must use same parameters to validate
3. **Rescue threshold:** With 7 seeds, threshold = 3-4 (not 2-3)

---

## ALGORITHM 1 CHECK

Paper provides Algorithm 1 (vectorized chaining). We implement a simpler version:

**Paper (Algorithm 1):**
- Sorts all anchors by position
- Merges anchors within C_b tolerance
- Vectorized with NumPy/PyTorch
- O(N·MK·log(MK)) complexity

**Our Implementation:**
- Groups anchors by expected position
- Merges within tolerance C
- Python loops (not vectorized)
- Similar logic, less optimized

**Status:** ⚠️ **FUNCTIONALLY EQUIVALENT** but not optimized

---

## VERDICT

### What's Correct ✅
- Core chaining algorithm (Equation 2)
- Index configuration (A.3)
- Training protocol (A.2)
- Architecture design
- Score calculation

### What's Wrong ❌
- **Seed numbers: 5→16 instead of 7→13**

### What's Unspecified ⚠️
- Tolerance C value
- Top-K chains
- Temperature τ

**ACTION REQUIRED:** Change seed numbers from 5→16 to 7→13 to match paper exactly.
