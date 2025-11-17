import argparse, time, json, subprocess, os, numpy as np, torch, faiss
from pathlib import Path
from tqdm import tqdm
from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode

def run_wfa(query, ref, wfa_bin, window):
    tmp_q = "/tmp/q.fa"
    tmp_r = "/tmp/r.fa"
    with open(tmp_q, "w") as f: f.write(">q\n" + query + "\n")
    with open(tmp_r, "w") as f: f.write(">r\n" + ref + "\n")
    cmd = [wfa_bin, "-a", "gap-affine-wfa", "-i", tmp_q, "-o", "/tmp/wfa_out.txt"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not os.path.exists("/tmp/wfa_out.txt"): return None
    with open("/tmp/wfa_out.txt") as f:
        for line in f:
            if line and not line.startswith("["):
                return line.strip()
    return None

def main(args):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model,_ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()

    idx = faiss.read_index(args.index)
    positions = np.load(args.positions)
    ref = ReferenceGenome(args.fasta)
    seq = ref.sequences[args.chrom]

    rng = np.random.RandomState(args.seed)
    reads = []
    true_pos = []
    for _ in range(args.num_reads):
        p = rng.randint(0, len(seq)-args.read_len)
        reads.append(seq[p:p+args.read_len])
        true_pos.append(p)

    stats = {"total": args.num_reads, "mapped":0, "unmapped":0, "no_candidate":0, "wfa_fail":0}

    for r,pos in tqdm(zip(reads,true_pos), total=len(reads), desc="reads"):
        seed = r[:model.input_len]
        t = torch.tensor(one_hot_encode(seed)).unsqueeze(0).to(device).float()
        with torch.no_grad():
            emb = model(t).cpu().numpy()
        _, I = idx.search(emb, args.k)
        cands = positions[I[0]]
        if len(cands)==0:
            stats["no_candidate"]+=1
            stats["unmapped"]+=1
            continue

        ref_center = cands[0]
        ref_slice = seq[max(0,ref_center-args.window): ref_center+args.window]
        cigar = run_wfa(r, ref_slice, args.wfa_bin, args.window)
        if cigar:
            stats["mapped"]+=1
        else:
            stats["wfa_fail"]+=1
            stats["unmapped"]+=1

    print(json.dumps({"stats":stats}, indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--checkpoint",required=True)
    p.add_argument("--index",required=True)
    p.add_argument("--positions",required=True)
    p.add_argument("--fasta",required=True)
    p.add_argument("--chrom",required=True)
    p.add_argument("--num-reads",type=int,default=100)
    p.add_argument("--read-len",type=int,default=2000)
    p.add_argument("--k",type=int,default=10)
    p.add_argument("--window",type=int,default=8192)
    p.add_argument("--seed",type=int,default=42)
    p.add_argument("--wfa-bin",default=os.environ.get("WFA_BIN","align_benchmark"))
    args=p.parse_args()
    main(args)
