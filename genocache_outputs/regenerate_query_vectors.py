import numpy as np, torch, json, os, sys
from Bio import SeqIO
import torch.nn as nn

# ---------- params (match training) ----------
SEED_LEN = 256
OFFSETS = [0,128,256,384,512,640]
BATCH = 128
OUT_QV = "genocache_outputs/query_vectors.npy"
OUT_IDS = "genocache_outputs/query_ids.json"
OUT_OFF = "genocache_outputs/query_offsets.npy"
CHECKPOINT = "nal_encoder_epoch1.pt"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
READS = "genocache_data/reads_chr22_synth_1kb_500.fa"

class TinyNAL(nn.Module):
    def __init__(self, out_dim=128):
        super().__init__()
        self.in_proj = nn.Conv1d(4, 64, kernel_size=7, padding=3)
        self.conv = nn.Sequential(
            nn.ReLU(),
            nn.Conv1d(64, 128, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.proj = nn.Linear(128, out_dim)
    def forward(self, x):
        x = x.permute(0,2,1)
        x = self.in_proj(x)
        x = self.conv(x)
        x = x.view(x.size(0), -1)
        x = self.proj(x)
        return nn.functional.normalize(x, dim=-1)

def seq_to_onehot(seq):
    arr = np.zeros((len(seq),4), dtype=np.float32)
    for i,ch in enumerate(seq):
        if ch=='A': arr[i,0]=1
        elif ch=='C': arr[i,1]=1
        elif ch=='G': arr[i,2]=1
        elif ch=='T': arr[i,3]=1
        else: arr[i,np.random.randint(4)] = 1
    return arr

model = TinyNAL(128)
model.load_state_dict(torch.load(CHECKPOINT, map_location="cpu"))
model.to(DEVICE).eval()

seeds = []
meta = []
for rec in SeqIO.parse(READS, "fasta"):
    s = str(rec.seq).upper()
    rid = rec.id
    L = len(s)
    for off in OFFSETS:
        seg = s[off:off+SEED_LEN]
        if len(seg) < SEED_LEN:
            seg += "A"*(SEED_LEN-len(seg))
        oh = seq_to_onehot(seg)
        seeds.append(oh)
        meta.append((rid, off, L))

print("Total query seeds:", len(seeds), "batch size:", BATCH, "device:", DEVICE)

vecs = []
for i in range(0, len(seeds), BATCH):
    batch = np.stack(seeds[i:i+BATCH], axis=0).astype("float32")
    batch = np.ascontiguousarray(batch).copy()  # ✅ FIX: force clean contiguous copy
    t = torch.from_numpy(batch).to(DEVICE)
    with torch.no_grad():
        v = model(t).cpu().numpy()
    vecs.append(v)

vecs = np.vstack(vecs)
norms = np.linalg.norm(vecs, axis=1, keepdims=True)
norms[norms==0]=1
vecs = vecs / norms

np.save(OUT_QV, vecs)
json.dump([m[0] for m in meta], open(OUT_IDS,"w"))
np.save(OUT_OFF, np.array([m[1] for m in meta], dtype=np.int32))

print("Saved query vectors:", OUT_QV, vecs.shape)
print("Saved query ids:", OUT_IDS)
print("Saved offsets :", OUT_OFF)
