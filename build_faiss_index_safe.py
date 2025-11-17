#!/usr/bin/env python3
import os, sys, math
import numpy as np
import faiss

REF_VECTORS = "ref_vectors.npy"
OUT_INDEX = "faiss_index.ivfpq"
OUT_INDEX_FLAT = "faiss_index_flat.ivf"
NPROBE = 8
PQ_M = 16

if not os.path.exists(REF_VECTORS):
    raise SystemExit("Missing ref_vectors.npy - run encode_ref_vectors.py first")

ref_vecs = np.load(REF_VECTORS).astype("float32")
N, DIM = ref_vecs.shape
print("Loaded ref_vectors:", REF_VECTORS, "shape:", ref_vecs.shape)

# small dataset => exact index
if N < 2048:
    print("N small; building exact IndexFlatIP (no training).")
    index = faiss.IndexFlatIP(DIM)
    index.add(ref_vecs)
    index.nprobe = min(N, NPROBE)
    faiss.write_index(index, OUT_INDEX_FLAT)
    print("Saved exact index:", OUT_INDEX_FLAT, "ntotal=", index.ntotal)
    sys.exit(0)

# otherwise: safe IVFPQ
nlist = max(32, int(math.sqrt(N)))
nlist = min(nlist, N)
print("Using nlist=", nlist, "nprobe=", NPROBE, "PQ_M=", PQ_M)

quantizer = faiss.IndexFlatIP(DIM)
index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)

# choose train_size >= nlist and ideally >= 5*nlist
train_size = min(max(nlist * 5, 10000), N)
if train_size < nlist:
    train_size = nlist
print("Training size:", train_size, "of", N)

if not index.is_trained:
    print("Training IVFPQ index on first", train_size, "vectors...")
    index.train(ref_vecs[:train_size].copy())
    print("Training complete.")

print("Adding", N, "vectors to index...")
index.add(ref_vecs)
print("Index built. ntotal =", index.ntotal)

index.nprobe = NPROBE

# try to move to GPU if faiss GPU available
try:
    res = faiss.StandardGpuResources()
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    # save CPU index as canonical file
    faiss.write_index(index, OUT_INDEX)
    print("Moved to GPU for runtime. Saved CPU index:", OUT_INDEX)
except Exception as e:
    print("GPU Faiss not available or failed, staying on CPU. Reason:", e)
    faiss.write_index(index, OUT_INDEX)
    print("Saved CPU index:", OUT_INDEX)

print("Done.")
