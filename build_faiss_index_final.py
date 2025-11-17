#!/usr/bin/env python3
import os, sys, math, numpy as np

REF_VECTORS = "ref_vectors.npy"
REF_POSITIONS = "ref_positions.npy"
OUT_INDEX_CPU = "faiss_index_cpu.ivf"
OUT_INDEX_GPU = "faiss_index_gpu.ivf"

MIN_NLIST = 32
NPROBE = 8
PQ_M = 16

if not os.path.exists(REF_VECTORS):
    print("Missing", REF_VECTORS); sys.exit(1)

ref_vecs = np.load(REF_VECTORS).astype("float32")
N, DIM = ref_vecs.shape
print("Loaded ref_vectors:", REF_VECTORS, "shape:", ref_vecs.shape)

if N < 2048:
    print("Small dataset; building exact IndexFlatIP (no training).")
    import faiss
    index = faiss.IndexFlatIP(DIM)
    index.add(ref_vecs)
    faiss.write_index(index, OUT_INDEX_CPU)
    print("Wrote exact IndexFlatIP ->", OUT_INDEX_CPU)
    sys.exit(0)

nlist = max(MIN_NLIST, int(math.sqrt(N)))
nlist = min(nlist, N)
print("nlist (clusters) chosen:", nlist)

if DIM % PQ_M != 0:
    # choose divisor <= PQ_M
    for m in range(PQ_M, 0, -1):
        if DIM % m == 0:
            PQ_M = m
            break
    else:
        PQ_M = 1
print("Using PQ_M =", PQ_M, " (DIM % PQ_M == 0 ? ", (DIM % PQ_M == 0), ")")

train_size = min(100_000, N)
train_size = max(train_size, 5 * nlist)
if train_size > N:
    train_size = N
print("train_size for clustering:", train_size)

try:
    import faiss
    quantizer = faiss.IndexFlatIP(DIM)
    index = faiss.IndexIVFPQ(quantizer, DIM, nlist, PQ_M, 8)
    index.metric_type = faiss.METRIC_INNER_PRODUCT
    if not index.is_trained:
        print("Training IVFPQ index...")
        rnd = np.random.RandomState(1234)
        perm = rnd.permutation(N)[:train_size]
        index.train(ref_vecs[perm])
    print("Adding vectors...")
    index.add(ref_vecs)
    index.nprobe = NPROBE
    faiss.write_index(index, OUT_INDEX_CPU)
    print("Saved CPU index:", OUT_INDEX_CPU)
except Exception as e:
    print("Faiss IVFPQ build failed — falling back to IndexFlatIP. Reason:", e)
    try:
        import faiss
        index = faiss.IndexFlatIP(DIM)
        index.add(ref_vecs)
        faiss.write_index(index, OUT_INDEX_CPU)
        print("Saved fallback exact index:", OUT_INDEX_CPU)
    except Exception as e2:
        print("Fallback exact index build also failed:", e2)
        raise

# Optional: try to move to GPU if GPU Faiss is available
try:
    import faiss
    if hasattr(faiss, "StandardGpuResources"):
        print("Attempting to move index to GPU.")
        res = faiss.StandardGpuResources()
        gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
        gpu_index.nprobe = NPROBE
        faiss.write_index(faiss.index_gpu_to_cpu(gpu_index), OUT_INDEX_GPU)
        print("Moved to GPU and saved CPU-copy:", OUT_INDEX_GPU)
    else:
        print("GPU Faiss bindings not present; staying CPU.")
except Exception as e:
    print("Optional GPU move failed (nonfatal). Reason:", e)

print("Done.")
