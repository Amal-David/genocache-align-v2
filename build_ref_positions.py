#!/usr/bin/env python3
import os
import re
import numpy as np
from Bio import SeqIO

FA = "chr22_1kb.fa"
OUT = "ref_positions.npy"
STRIDE = None   # set to int to override, else inferred if headers contain coordinates

if not os.path.exists(FA):
    raise SystemExit(f"Missing {FA} in cwd {os.getcwd()}")

positions = []
for rec in SeqIO.parse(FA, "fasta"):
    h = rec.description
    # try to parse coordinates like chr22:12345-13344
    m = re.search(r'[:](\d+)-(\d+)', h)
    if m:
        start = int(m.group(1))
    else:
        positions.append(None)
        continue
    positions.append(start)

if STRIDE is not None:
    # override: compute positions as 0, stride, 2*stride ...
    N = len(list(SeqIO.parse(FA, "fasta")))
    positions = [i*STRIDE for i in range(N)]

# if any header failed parse, fallback to synthetic positions with inferred stride = len(seq) - overlap unknown
if any(p is None for p in positions):
    print("Header parsing failed for some records; building synthetic start positions with stride=32")
    seqs = list(SeqIO.parse(FA, "fasta"))
    N = len(seqs)
    stride = 32
    positions = [i*stride for i in range(N)]

positions = np.array(positions, dtype=np.int64)
np.save(OUT, positions)
print("Wrote", OUT, "len=", len(positions))
