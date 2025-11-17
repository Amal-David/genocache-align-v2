#!/usr/bin/env python3
"""
Benchmark encoding and search speed
"""
import time, random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import faiss

SEED_LEN = 256
EMB_DIM = 256
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

class ImprovedNAL(nn.Module):
    def __init__(self, out_dim=EMB_DIM):
        super().__init__()
        self.in_proj = nn.Conv1d(4, 128, kernel_size=9, padding=4)
        self.bn1 = nn.BatchNorm1d(128)
        self.conv1 = nn.Conv1d(128, 256, kernel_size=7, padding=3)
        self.bn2 = nn.BatchNorm1d(256)
        self.conv2 = nn.Conv1d(256, 256, kernel_size=7, padding=3)
        self.bn3 = nn.BatchNorm1d(256)
        self.conv3 = nn.Conv1d(256, 512, kernel_size=5, padding=2)
        self.bn4 = nn.BatchNorm1d(512)
        self.conv4 = nn.Conv1d(512, 512, kernel_size=5, padding=2)
        self.bn5 = nn.BatchNorm1d(512)
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.proj = nn.Sequential(nn.Linear(512, out_dim), nn.LayerNorm(out_dim))
    
    def forward(self, x):
        x = x.permute(0,2,1)
        x = F.relu(self.bn1(self.in_proj(x)))
        identity = x
        x = F.relu(self.bn2(self.conv1(x)))
        x = self.bn3(self.conv2(x))
        if identity.shape[1] != x.shape[1]:
            identity = F.conv1d(identity, torch.eye(256, 128, device=x.device).unsqueeze(2), padding=0)
        x = F.relu(x + identity)
        x = F.relu(self.bn4(self.conv3(x)))
        identity = x
        x = self.bn5(self.conv4(x))
        x = F.relu(x + identity)
        x = self.pool(x).view(x.shape[0], -1)
        x = self.proj(x)
        return F.normalize(x, dim=-1)

def seq_to_onehot(seq):
    arr = np.zeros((len(seq),4), dtype=np.float32)
    for i,ch in enumerate(seq):
        ch = ch.upper()
        if ch == "A": arr[i,0]=1
        elif ch=="C": arr[i,1]=1
        elif ch=="G": arr[i,2]=1
        elif ch=="T": arr[i,3]=1
        elif ch=="N": arr[i, random.randrange(4)] = 0.25
        else: arr[i, random.randrange(4)] = 1.0
    return arr

def generate_random_sequences(n, length=SEED_LEN):
    """Generate random DNA sequences"""
    bases = ['A', 'C', 'G', 'T']
    seqs = []
    for _ in range(n):
        seq = ''.join(random.choices(bases, k=length))
        seqs.append(seq)
    return seqs

def benchmark_encoding(model, batch_sizes=[1, 10, 100, 256, 512, 1024], n_trials=10):
    """Benchmark encoding speed"""
    print("="*70)
    print("Encoding Speed Benchmark")
    print("="*70)
    
    print(f"\nDevice: {DEVICE}")
    print(f"Sequence length: {SEED_LEN}bp")
    print(f"Embedding dimension: {EMB_DIM}")
    print(f"Trials per batch size: {n_trials}")
    
    print(f"\n{'Batch Size':>12} {'Time (ms)':>12} {'Seqs/sec':>12} {'Throughput':>15}")
    print("-" * 53)
    
    results = []
    
    for batch_size in batch_sizes:
        # Generate test data
        seqs = generate_random_sequences(batch_size)
        arr = np.stack([seq_to_onehot(s) for s in seqs])
        tensor = torch.from_numpy(arr).to(DEVICE)
        
        # Warmup
        with torch.no_grad():
            for _ in range(5):
                _ = model(tensor)
        
        # Benchmark
        times = []
        for _ in range(n_trials):
            torch.cuda.synchronize() if torch.cuda.is_available() else None
            start = time.time()
            
            with torch.no_grad():
                _ = model(tensor)
            
            torch.cuda.synchronize() if torch.cuda.is_available() else None
            elapsed = time.time() - start
            times.append(elapsed)
        
        avg_time = np.mean(times) * 1000  # ms
        seqs_per_sec = batch_size / (np.mean(times))
        mbp_per_sec = (batch_size * SEED_LEN) / (np.mean(times)) / 1e6
        
        print(f"{batch_size:>12} {avg_time:>11.2f} {seqs_per_sec:>11.0f} {mbp_per_sec:>10.2f} Mbp/s")
        
        results.append({
            'batch_size': batch_size,
            'time_ms': avg_time,
            'seqs_per_sec': seqs_per_sec,
            'mbp_per_sec': mbp_per_sec
        })
    
    # Find optimal batch size
    best = max(results, key=lambda x: x['seqs_per_sec'])
    print(f"\n✓ Optimal batch size: {best['batch_size']} ({best['seqs_per_sec']:.0f} seqs/sec)")
    
    return results

def benchmark_search(index, query_sizes=[10, 100, 1000, 10000], k_values=[1, 10, 20], n_trials=10):
    """Benchmark FAISS search speed"""
    print("\n" + "="*70)
    print("FAISS Search Speed Benchmark")
    print("="*70)
    
    print(f"\nIndex size: {index.ntotal:,} vectors")
    print(f"Index type: {type(index).__name__}")
    if hasattr(index, 'nprobe'):
        print(f"nprobe: {index.nprobe}")
    
    results = []
    
    for n_queries in query_sizes:
        # Generate random query vectors
        query_vecs = np.random.randn(n_queries, EMB_DIM).astype('float32')
        query_vecs = query_vecs / np.linalg.norm(query_vecs, axis=1, keepdims=True)
        
        for k in k_values:
            # Warmup
            for _ in range(5):
                _ = index.search(query_vecs, k)
            
            # Benchmark
            times = []
            for _ in range(n_trials):
                start = time.time()
                D, I = index.search(query_vecs, k)
                elapsed = time.time() - start
                times.append(elapsed)
            
            avg_time = np.mean(times) * 1000  # ms
            queries_per_sec = n_queries / np.mean(times)
            us_per_query = np.mean(times) * 1e6 / n_queries
            
            results.append({
                'n_queries': n_queries,
                'k': k,
                'time_ms': avg_time,
                'queries_per_sec': queries_per_sec,
                'us_per_query': us_per_query
            })
    
    print(f"\n{'Queries':>10} {'K':>5} {'Time (ms)':>12} {'Queries/sec':>13} {'μs/query':>12}")
    print("-" * 54)
    
    for r in results:
        print(f"{r['n_queries']:>10} {r['k']:>5} {r['time_ms']:>11.2f} "
              f"{r['queries_per_sec']:>12.0f} {r['us_per_query']:>11.2f}")
    
    # Best performance
    best = max(results, key=lambda x: x['queries_per_sec'])
    print(f"\n✓ Best performance: {best['queries_per_sec']:.0f} queries/sec "
          f"(n={best['n_queries']}, k={best['k']})")
    
    return results

def benchmark_end_to_end(model, index, n_reads=1000, batch_size=256, k=20):
    """Benchmark complete alignment pipeline"""
    print("\n" + "="*70)
    print("End-to-End Pipeline Benchmark")
    print("="*70)
    
    print(f"\nConfiguration:")
    print(f"  Reads: {n_reads:,}")
    print(f"  Batch size: {batch_size}")
    print(f"  Top-K: {k}")
    
    # Generate sequences
    print(f"\n[1/3] Generating {n_reads:,} random sequences...")
    start = time.time()
    seqs = generate_random_sequences(n_reads)
    gen_time = time.time() - start
    print(f"  ✓ Generated in {gen_time:.2f}s")
    
    # Encode
    print(f"\n[2/3] Encoding...")
    all_vecs = []
    encode_start = time.time()
    
    with torch.no_grad():
        for i in range(0, n_reads, batch_size):
            batch = seqs[i:i+batch_size]
            arr = np.stack([seq_to_onehot(s) for s in batch])
            tensor = torch.from_numpy(arr).to(DEVICE)
            vecs = model(tensor).cpu().numpy()
            all_vecs.append(vecs)
    
    query_vecs = np.vstack(all_vecs).astype('float32')
    encode_time = time.time() - encode_start
    print(f"  ✓ Encoded in {encode_time:.2f}s ({n_reads/encode_time:.0f} reads/sec)")
    
    # Search
    print(f"\n[3/3] Searching...")
    search_start = time.time()
    D, I = index.search(query_vecs, k)
    search_time = time.time() - search_start
    print(f"  ✓ Searched in {search_time:.2f}s ({n_reads/search_time:.0f} queries/sec)")
    
    # Summary
    total_time = encode_time + search_time
    print(f"\n{'='*70}")
    print("Summary:")
    print(f"  Encoding:   {encode_time:>8.2f}s ({encode_time/total_time*100:5.1f}%)")
    print(f"  Searching:  {search_time:>8.2f}s ({search_time/total_time*100:5.1f}%)")
    print(f"  Total:      {total_time:>8.2f}s")
    print(f"  Throughput: {n_reads/total_time:>8.0f} reads/sec")
    print(f"  Latency:    {total_time/n_reads*1000:>8.2f} ms/read")
    print("="*70)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark GenoCache performance')
    parser.add_argument('--model', default='nal_encoder_best.pt', help='Model checkpoint')
    parser.add_argument('--index', default='faiss_index_improved.idx', help='FAISS index')
    parser.add_argument('--quick', action='store_true', help='Quick benchmark (fewer trials)')
    args = parser.parse_args()
    
    print("GenoCache Performance Benchmark")
    print("="*70)
    
    # Load model
    print(f"\nLoading model: {args.model}")
    model = ImprovedNAL().to(DEVICE)
    model.load_state_dict(torch.load(args.model, map_location=DEVICE))
    model.eval()
    print(f"  ✓ Loaded on {DEVICE}")
    
    # Load index
    print(f"\nLoading index: {args.index}")
    index = faiss.read_index(args.index)
    print(f"  ✓ Loaded {index.ntotal:,} vectors")
    
    # Run benchmarks
    n_trials = 3 if args.quick else 10
    
    # Encoding benchmark
    batch_sizes = [1, 10, 100, 256, 512] if args.quick else [1, 10, 100, 256, 512, 1024]
    encoding_results = benchmark_encoding(model, batch_sizes, n_trials)
    
    # Search benchmark
    query_sizes = [100, 1000] if args.quick else [10, 100, 1000, 10000]
    k_values = [1, 10, 20]
    search_results = benchmark_search(index, query_sizes, k_values, n_trials)
    
    # End-to-end benchmark
    n_reads = 500 if args.quick else 1000
    benchmark_end_to_end(model, index, n_reads)
    
    print("\n✓ Benchmark complete!")

if __name__ == '__main__':
    main()
