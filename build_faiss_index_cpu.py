#!/usr/bin/env python3
import os, math
import numpy as np
import faiss

# —— config (tweak if needed) ——
REF_VECTORS = "ref_vectors.npy"     # produced earlier in ~/genocache
QUERY_VECTORS = "query_vectors.npy" # optional; if missing we'll query the refs
TOPK = 10

# —— load ref vectors ——
if not os.path.exists(REF_VECTORS):
    raise SystemExit(f"Missing {REF_VECTORS} in cwd {os.getcwd()}")

ref_vecs = np.load(REF_VECTORS).astype("float32")
N, DIM = ref_vecs.shape
print("Loaded ref_vectors:", REF_VECTORS, "shape:", ref_vecs.shape)

# normalize (if embeddings are intended to be cosine / inner product)
def maybe_normalize(x):
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return x / norms

ref_vecs = maybe_normalize(ref_vecs)

# sensible defaults
nlist = max(1, int(math.sqrt(N)))
nprobe = 8
PQ_M = 16   # PQ16x8 -> M=16

print("Using nlist=", nlist, "nprobe=", nprobe, "PQ_M=", PQ_M)

# Build CPU IVFPQ index
quantizer = faiss.IndexFlatIP(DIM)  # inner product; cosine if normalized
index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)  # 8 = bits per subvector code

# train if necessary
if not index.is_trained:
    print("Training IVFPQ index on reference vectors...")
    index.train(ref_vecs)
    print("Training done.")

print("Adding vectors to index...")
index.add(ref_vecs)
print("Index built. ntotal =", index.ntotal)

index.nprobe = nprobe

# quick test: search the first ref vector
q = ref_vecs[:1]
D, I = index.search(q, TOPK)
print("Top{} for first ref vector:".format(TOPK))
print("Indices:", I)
print("Scores :", D)

# optional: run batch queries if query_vectors exists
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
