#!/usr/bin/env python3
import os, math, sys
import numpy as np
try:
    import faiss
except Exception as e:
    raise SystemExit("import faiss failed: " + str(e))

REF_VECTORS = "ref_vectors.npy"     # must be in cwd
QUERY_VECTORS = "query_vectors.npy" # optional
TOPK = 10

if not os.path.exists(REF_VECTORS):
    raise SystemExit(f"Missing {REF_VECTORS} in cwd {os.getcwd()}")

ref_vecs = np.load(REF_VECTORS).astype("float32")
N, DIM = ref_vecs.shape
print("Loaded ref_vectors:", REF_VECTORS, "shape:", ref_vecs.shape)

# normalize (cosine-like behaviour for inner product)
def maybe_normalize(x):
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return x / norms

ref_vecs = maybe_normalize(ref_vecs)

# Heuristic: small datasets -> use exact flat index (no training)
SMALL_THRESHOLD = 2048   # below this, use brute-force
if N < SMALL_THRESHOLD:
    print(f"N={N} < {SMALL_THRESHOLD}: building exact IndexFlatIP (no training).")
    index = faiss.IndexFlatIP(DIM)
    index.add(ref_vecs)
    print("IndexFlatIP built. ntotal =", index.ntotal)
else:
    # safe nlist: at most N (and typically sqrt(N)), but cap it reasonably
    nlist = max(1, int(math.sqrt(N)))
    # make sure nlist < N
    if nlist >= N:
        nlist = max(1, N // 2)
    nlist = min(nlist, N // 2) if N >= 4 else 1
    nprobe = 8
    PQ_M = 16
    print("Attempting IVFPQ with nlist=", nlist, "nprobe=", nprobe, "PQ_M=", PQ_M)
    quantizer = faiss.IndexFlatIP(DIM)
    index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)
    try:
        if not index.is_trained:
            print("Training IVFPQ index...")
            index.train(ref_vecs)
            print("Training done.")
        index.add(ref_vecs)
        index.nprobe = nprobe
        print("IVFPQ index built. ntotal =", index.ntotal)
    except Exception as e:
        print("IVFPQ train/add failed:", e)
        print("Falling back to IndexFlatIP (exact index).")
        index = faiss.IndexFlatIP(DIM)
        index.add(ref_vecs)
        print("IndexFlatIP built. ntotal =", index.ntotal)

# quick search test: first vector
q = ref_vecs[:1]
D, I = index.search(q, TOPK)
print("Top{} for first ref vector:".format(TOPK))
print("Indices:", I)
print("Scores :", D)

if os.path.exists(QUERY_VECTORS):
    qvs = np.load(QUERY_VECTORS).astype("float32")
    qvs = maybe_normalize(qvs)
    print("Loaded query vectors:", QUERY_VECTORS, "shape:", qvs.shape)
    Dq, Iq = index.search(qvs, TOPK)
    print("Sample query results (first 5):")
    for i in range(min(5, qvs.shape[0])):
        print("q{}: top idxs {} scores {}".format(i, Iq[i].tolist(), [float(x) for x in Dq[i].tolist()]))
else:
    print(f"No {QUERY_VECTORS} found; skipped batch query test.")

print("Done.")
