#!/usr/bin/env python3
import argparse
import time
import numpy as np
import faiss
import torch
from pathlib import Path
from Bio import SeqIO
import sys
import subprocess
import re

sys.path.insert(0, str(Path(__file__).parent))
from improved_cnn import ImprovedCNN
from genomic_utils import one_hot_encode

def load_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences

def cluster_positions(positions, window=100):
    if len(positions) == 0:
        return []
    positions = sorted(positions)
    clusters = []
    current_cluster = [positions[0]]
    for pos in positions[1:]:
        if pos - current_cluster[-1] <= window:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)
    return clusters

def wfa_to_sam_cigar(wfa_cigar):
    """Convert WFA CIGAR (MMMMXMMM) to SAM CIGAR (4M1X3M)"""
    if not wfa_cigar:
        return "*"
    
    sam_ops = []
    current_op = wfa_cigar[0]
    count = 1
    
    for op in wfa_cigar[1:]:
        if op == current_op:
            count += 1
        else:
            sam_ops.append(f"{count}{current_op}")
            current_op = op
            count = 1
    sam_ops.append(f"{count}{current_op}")
    
    return ''.join(sam_ops)

def align_with_wfa(read_seq, ref_seq, wfa_bin):
    """Use WFA2 C++ tool to get alignment"""
    try:
        result = subprocess.run(
            [wfa_bin, read_seq, ref_seq],
            capture_output=True, text=True, timeout=5
        )
        
        if result.returncode == 0:
            # Parse output: "score\tCIGAR"
            parts = result.stdout.strip().split('\t')
            if len(parts) == 2:
                score = int(parts[0])
                wfa_cigar = parts[1]
                sam_cigar = wfa_to_sam_cigar(wfa_cigar)
                return sam_cigar, score, True
        
        return "*", 0, False
    except Exception as e:
        return "*", 0, False

def extract_seeds_batch(reads, seed_len=512, num_seeds=5):
    all_seeds = []
    read_seed_map = []
    
    for read_id, read_seq in enumerate(reads):
        read_seeds = []
        for i in range(num_seeds):
            offset = i * ((len(read_seq) - seed_len) // 4)
            seed = read_seq[offset:offset+seed_len]
            all_seeds.append((seed, offset))
            read_seeds.append(len(all_seeds) - 1)
        read_seed_map.append(read_seeds)
    
    return all_seeds, read_seed_map

def encode_seeds_batch(model, seeds, batch_size=512, device='cuda'):
    embeddings = []
    
    for i in range(0, len(seeds), batch_size):
        batch_seeds = seeds[i:i+batch_size]
        batch_tensor = torch.zeros((len(batch_seeds), 4, 512), dtype=torch.float32)
        
        for j, (seed, _) in enumerate(batch_seeds):
            try:
                encoded = one_hot_encode(seed)
                batch_tensor[j] = torch.tensor(encoded, dtype=torch.float32)
            except:
                pass
        
        batch_tensor = batch_tensor.to(device)
        
        with torch.no_grad():
            batch_emb = model(batch_tensor).cpu().numpy()
        
        embeddings.extend(batch_emb)
    
    return np.array(embeddings, dtype=np.float32)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--num-reads", type=int, default=100)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--output", default="aligned.sam")
    p.add_argument("--window-size", type=int, default=500, help="Reference window padding")
    args = p.parse_args()
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    wfa_bin = str(Path(__file__).parent.parent / "wfa_align")
    
    print(f"🚀 Full GenoCache-Align Pipeline")
    print(f"   Device: {device}")
    print(f"   WFA2: {wfa_bin}")
    
    # Load everything
    start = time.time()
    model, _ = ImprovedCNN.load_checkpoint(args.checkpoint, device=device)
    model.eval()
    
    index = faiss.read_index(args.index)
    ref_positions = np.load(args.positions)
    sequences = load_fasta(args.fasta)
    ref_seq = sequences[args.chrom]
    
    load_time = time.time() - start
    print(f"✅ Loaded in {load_time:.1f}s")
    
    # Generate reads
    print(f"\n📝 Generating {args.num_reads} reads...")
    np.random.seed(42)
    reads = []
    true_positions = []
    
    max_start = len(ref_seq) - 2000 - 10
    candidate_starts = np.random.randint(10, max_start, args.num_reads)
    
    for pos in candidate_starts:
        read = ref_seq[pos:pos+2000]
        if len(read) == 2000 and 'N' not in read:
            reads.append(read)
            true_positions.append(pos)
    
    print(f"   Valid reads: {len(reads)}")
    
    # Seed-based mapping
    print(f"\n🧬 Seed mapping & encoding...")
    map_start = time.time()
    all_seeds, read_seed_map = extract_seeds_batch(reads)
    all_embeddings = encode_seeds_batch(model, all_seeds, args.batch_size, device)
    D, I = index.search(all_embeddings, k=10)
    map_time = time.time() - map_start
    
    # Vote for positions
    print(f"\n🗳️  Voting for positions...")
    vote_start = time.time()
    mapped_reads = []
    
    for read_idx, seed_indices in enumerate(read_seed_map):
        all_positions = []
        for seed_idx in seed_indices:
            indices = I[seed_idx]
            positions = ref_positions[indices]
            seed_offset = all_seeds[seed_idx][1]
            adjusted = positions - seed_offset
            all_positions.extend(adjusted)
        
        clusters = cluster_positions(all_positions, window=100)
        best_cluster = max(clusters, key=len)
        predicted_pos = int(np.median(best_cluster))
        
        mapped_reads.append({
            'read_id': f"read_{read_idx}",
            'seq': reads[read_idx],
            'pos': predicted_pos,
            'true_pos': true_positions[read_idx]
        })
    
    vote_time = time.time() - vote_start
    
    # WFA2 alignment
    print(f"\n🔬 WFA2 fine alignment...")
    align_start = time.time()
    sam_records = []
    correct = 0
    
    for i, read_info in enumerate(mapped_reads):
        pos = read_info['pos']
        
        # Extract reference window
        ref_start = max(0, pos - args.window_size)
        ref_end = min(len(ref_seq), pos + len(read_info['seq']) + args.window_size)
        ref_window = ref_seq[ref_start:ref_end]
        
        # Align with WFA2
        cigar, score, success = align_with_wfa(read_info['seq'], ref_window, wfa_bin)
        
        # Check accuracy (within 1kb)
        if abs(pos - read_info['true_pos']) <= 1024:
            correct += 1
        
        # Create SAM record
        flag = 0 if success else 4
        mapq = 60 if success else 0
        sam_records.append(
            f"{read_info['read_id']}\t{flag}\t{args.chrom}\t{pos+1}\t{mapq}\t{cigar}\t*\t0\t0\t{read_info['seq']}\t*\tAS:i:{score}"
        )
        
        if (i+1) % 50 == 0:
            print(f"   Aligned {i+1}/{len(mapped_reads)} reads...")
    
    align_time = time.time() - align_start
    
    # Write SAM
    print(f"\n💾 Writing SAM to {args.output}...")
    with open(args.output, 'w') as f:
        f.write(f"@HD\tVN:1.6\tSO:unsorted\n")
        f.write(f"@SQ\tSN:{args.chrom}\tLN:{len(ref_seq)}\n")
        f.write(f"@PG\tID:genocache\tPN:GenoCache-Align\tVN:0.1\n")
        for record in sam_records:
            f.write(record + "\n")
    
    # Summary
    total_time = time.time() - start
    accuracy = 100 * correct / len(reads)
    throughput = len(reads) / total_time
    
    print(f"\n📊 Results:")
    print(f"   Reads: {len(reads)}")
    print(f"   Accuracy: {accuracy:.1f}% ({correct}/{len(reads)})")
    print(f"   Throughput: {throughput:.1f} reads/sec")
    print(f"\n⏱️  Timing:")
    print(f"   Load: {load_time:.2f}s")
    print(f"   Seed mapping: {map_time:.2f}s ({100*map_time/total_time:.0f}%)")
    print(f"   Voting: {vote_time:.2f}s ({100*vote_time/total_time:.0f}%)")
    print(f"   WFA alignment: {align_time:.2f}s ({100*align_time/total_time:.0f}%)")
    print(f"   Total: {total_time:.2f}s")
    print(f"\n✅ SAM file: {args.output}")

if __name__ == "__main__":
    main()
