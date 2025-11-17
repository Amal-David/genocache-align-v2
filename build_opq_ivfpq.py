import numpy as np, faiss
REF="ref_vectors.npy"
IDX="faiss_opq_ivfpq_M8_nprobe64.idx"
ref = np.load(REF).astype("float32")
N,D = ref.shape
nlist = max(32, int(N**0.5))
train_size = min(100000, N)
PQ_M = 8
print("N,D",N,D,"nlist",nlist,"train",train_size,"PQ_M",PQ_M)
opq = faiss.OPQMatrix(D, PQ_M)
opq.train(ref[:train_size])
ref_opq = opq.apply_py(ref)
q = faiss.IndexFlatIP(D)
idx = faiss.IndexIVFPQ(q, D, nlist, PQ_M, 8)
if not idx.is_trained:
    idx.train(ref_opq[:train_size])
idx.add(ref_opq)
idx.nprobe = 64
faiss.write_index(idx, IDX)
faiss.write_index(faiss.IndexPreTransform(opq, faiss.IndexFlatIP(D)), "opq_matrix.idx")
print("WROTE", IDX, "and opq_matrix.idx")
