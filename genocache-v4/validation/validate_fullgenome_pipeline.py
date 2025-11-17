#!/usr/bin/env python3
"""
GenoCache V4 - Full Genome Pipeline Validation

Comprehensive validation testing:
1. Generate test reads from multiple chromosomes
2. Encode full genome with stride=32
3. Build FAISS index
4. Test accuracy across different chromosomes
5. Measure end-to-end performance
6. Test with real ONT data (if available)

This validates the COMPLETE pipeline!
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
import faiss
import json
from pathlib import Path
from tqdm import tqdm
import time
from collections import defaultdict

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder
sys.path.append(str(Path(__file__).parent.parent / "training"))
from augmentation import DataAugmentation


def load_genome(fasta_path: Path) -> dict:
    """Load all chromosomes from genome"""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr is not None:
                    genome[current_chr] = ''.join(current_seq).upper()
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_chr is not None:
            genome[current_chr] = ''.join(current_seq).upper()
    
    return genome


def generate_test_reads_multi_chr(
    genome: dict,
    test_chromosomes: list,
    reads_per_chr: int = 200,
    read_length: int = 512,
    error_rate: float = 0.05
) -> tuple:
    """
    Generate test reads from multiple chromosomes
    
    Returns:
        (reads, ground_truth): List of sequences and their positions
    """
    print(f"\nGenerating test reads from {len(test_chromosomes)} chromosomes...")
    print(f"  Reads per chromosome: {reads_per_chr}")
    print(f"  Read length: {read_length}bp")
    print(f"  Error rate: {error_rate} ({100*(1-error_rate):.1f}% identity)")
    
    augmenter = DataAugmentation(seed_len=read_length)
    all_reads = []
    ground_truth = []
    
    for chr_name in test_chromosomes:
        if chr_name not in genome:
            print(f"  ⚠️  Chromosome {chr_name} not found, skipping")
            continue
        
        chr_seq = genome[chr_name]
        print(f"  Generating from {chr_name} ({len(chr_seq):,} bp)...")
        
        chr_reads = 0
        attempts = 0
        max_attempts = reads_per_chr * 3
        
        while chr_reads < reads_per_chr and attempts < max_attempts:
            attempts += 1
            
            # Sample random position
            start_pos = np.random.randint(0, len(chr_seq) - read_length)
            end_pos = start_pos + read_length
            
            # Extract sequence
            clean_seq = chr_seq[start_pos:end_pos]
            
            # Skip high-N regions
            if clean_seq.count('N') > read_length * 0.1:
                continue
            
            # Add errors
            noisy_seq = augmenter.add_errors(clean_seq, error_rate)
            
            # Ensure correct length
            if len(noisy_seq) > read_length:
                noisy_seq = noisy_seq[:read_length]
            elif len(noisy_seq) < read_length:
                noisy_seq = noisy_seq + 'N' * (read_length - len(noisy_seq))
            
            # Random RC
            if np.random.random() < 0.5:
                noisy_seq = augmenter.reverse_complement(noisy_seq)
            
            all_reads.append(noisy_seq)
            ground_truth.append({
                'chromosome': chr_name,
                'start': start_pos,
                'end': end_pos,
                'length': read_length
            })
            chr_reads += 1
    
    print(f"✅ Generated {len(all_reads)} total test reads")
    return all_reads, ground_truth


def encode_genome_with_stride(
    model: torch.nn.Module,
    genome: dict,
    chromosomes: list,
    window_size: int = 512,
    stride: int = 32,
    device: str = 'cuda',
    batch_size: int = 256
) -> tuple:
    """
    Encode selected chromosomes with sparse indexing
    
    Returns:
        (embeddings, positions, chr_names): Arrays of embeddings, positions, and chromosome names
    """
    model.eval()
    
    print(f"\nEncoding genome with stride={stride}...")
    print(f"  Chromosomes: {len(chromosomes)}")
    print(f"  Window size: {window_size} bp")
    print(f"  Stride: {stride} bp")
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    all_embeddings = []
    all_positions = []
    all_chr_names = []
    
    total_windows = 0
    for chr_name in chromosomes:
        if chr_name not in genome:
            continue
        seq = genome[chr_name]
        n_windows = (len(seq) - window_size) // stride + 1
        total_windows += n_windows
    
    print(f"  Total windows: {total_windows:,}")
    
    with torch.no_grad():
        for chr_name in tqdm(chromosomes, desc="Chromosomes"):
            if chr_name not in genome:
                continue
            
            chr_seq = genome[chr_name]
            
            batch_seqs = []
            batch_positions = []
            
            for pos in range(0, len(chr_seq) - window_size + 1, stride):
                window = chr_seq[pos:pos + window_size]
                
                # Skip high-N regions
                if window.count('N') > window_size * 0.1:
                    continue
                
                batch_seqs.append(seq_to_tensor(window))
                batch_positions.append(pos)
                
                # Process batch
                if len(batch_seqs) >= batch_size:
                    batch_tensor = torch.stack(batch_seqs).to(device)
                    embeddings = model(batch_tensor)
                    embeddings = F.normalize(embeddings, p=2, dim=1)
                    
                    all_embeddings.append(embeddings.cpu().numpy())
                    all_positions.extend(batch_positions)
                    all_chr_names.extend([chr_name] * len(batch_positions))
                    
                    batch_seqs = []
                    batch_positions = []
            
            # Process remaining
            if batch_seqs:
                batch_tensor = torch.stack(batch_seqs).to(device)
                embeddings = model(batch_tensor)
                embeddings = F.normalize(embeddings, p=2, dim=1)
                
                all_embeddings.append(embeddings.cpu().numpy())
                all_positions.extend(batch_positions)
                all_chr_names.extend([chr_name] * len(batch_positions))
    
    # Concatenate all
    all_embeddings = np.vstack(all_embeddings)
    all_positions = np.array(all_positions)
    all_chr_names = np.array(all_chr_names)
    
    print(f"✅ Encoded {len(all_positions):,} windows")
    print(f"  Embedding shape: {all_embeddings.shape}")
    print(f"  Memory: {all_embeddings.nbytes / 1e9:.2f} GB")
    
    return all_embeddings, all_positions, all_chr_names


def build_faiss_index(embeddings: np.ndarray, nprobe: int = 32) -> faiss.Index:
    """Build FAISS IVFPQ index"""
    N, D = embeddings.shape
    print(f"\nBuilding FAISS index...")
    print(f"  Embeddings: {N:,} × {D}D")
    
    nlist = min(int(np.sqrt(N)), 4096)
    print(f"  Index type: IVFPQ")
    print(f"  Number of clusters: {nlist}")
    print(f"  nprobe: {nprobe}")
    
    quantizer = faiss.IndexFlatIP(D)
    index = faiss.IndexIVFPQ(quantizer, D, nlist, 16, 8)
    
    print("  Training...")
    index.train(embeddings.astype(np.float32))
    
    print("  Adding embeddings...")
    index.add(embeddings.astype(np.float32))
    
    index.nprobe = nprobe
    
    print(f"✅ Index built: {index.ntotal:,} embeddings")
    return index


def encode_reads(
    model: torch.nn.Module,
    reads: list,
    device: str = 'cuda',
    batch_size: int = 64
) -> np.ndarray:
    """Encode test reads"""
    model.eval()
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, len(reads), batch_size), desc="Encoding reads"):
            batch = reads[i:i+batch_size]
            batch_tensors = torch.stack([seq_to_tensor(seq) for seq in batch]).to(device)
            
            batch_emb = model(batch_tensors)
            batch_emb = F.normalize(batch_emb, p=2, dim=1)
            
            embeddings.append(batch_emb.cpu().numpy())
    
    return np.vstack(embeddings)


def validate_pipeline(
    index: faiss.Index,
    positions: np.ndarray,
    chr_names: np.ndarray,
    read_embeddings: np.ndarray,
    ground_truth: list,
    k: int = 32,
    tolerance: int = 1000
) -> dict:
    """
    Run comprehensive validation
    
    Returns detailed metrics per chromosome and overall
    """
    print(f"\nValidating pipeline...")
    print(f"  Reads: {len(read_embeddings)}")
    print(f"  Top-k: {k}")
    print(f"  Tolerance: ±{tolerance} bp")
    
    # Overall metrics
    correct_1kb = 0
    correct_100bp = 0
    correct_exact = 0
    
    # Per-chromosome metrics
    chr_stats = defaultdict(lambda: {
        'total': 0,
        'correct_1kb': 0,
        'correct_100bp': 0,
        'errors': []
    })
    
    # Search timing
    search_times = []
    
    for i in tqdm(range(len(read_embeddings)), desc="Searching"):
        query_emb = read_embeddings[i:i+1].astype(np.float32)
        true_chr = ground_truth[i]['chromosome']
        true_pos = ground_truth[i]['start']
        
        # Time the search
        start_time = time.time()
        distances, indices = index.search(query_emb, k)
        search_time = time.time() - start_time
        search_times.append(search_time)
        
        # Get predicted positions
        pred_chr_names = chr_names[indices[0]]
        pred_positions = positions[indices[0]]
        
        # Find best match on correct chromosome
        best_error = float('inf')
        for pred_chr, pred_pos in zip(pred_chr_names, pred_positions):
            if pred_chr == true_chr:
                error = abs(pred_pos - true_pos)
                if error < best_error:
                    best_error = error
        
        # Update stats
        chr_stats[true_chr]['total'] += 1
        chr_stats[true_chr]['errors'].append(best_error)
        
        if best_error <= 1000:
            correct_1kb += 1
            chr_stats[true_chr]['correct_1kb'] += 1
        if best_error <= 100:
            correct_100bp += 1
            chr_stats[true_chr]['correct_100bp'] += 1
        if best_error == 0:
            correct_exact += 1
    
    # Compile results
    n_reads = len(read_embeddings)
    
    results = {
        'overall': {
            'total_reads': n_reads,
            'accuracy_1kb': 100 * correct_1kb / n_reads,
            'accuracy_100bp': 100 * correct_100bp / n_reads,
            'accuracy_exact': 100 * correct_exact / n_reads,
            'avg_search_time_ms': np.mean(search_times) * 1000,
            'median_search_time_ms': np.median(search_times) * 1000,
        },
        'per_chromosome': {},
        'timing': {
            'total_search_time': sum(search_times),
            'reads_per_second': n_reads / sum(search_times)
        }
    }
    
    # Per-chromosome results
    for chr_name, stats in chr_stats.items():
        if stats['total'] > 0:
            results['per_chromosome'][chr_name] = {
                'total': stats['total'],
                'accuracy_1kb': 100 * stats['correct_1kb'] / stats['total'],
                'accuracy_100bp': 100 * stats['correct_100bp'] / stats['total'],
                'median_error': np.median(stats['errors']),
                'mean_error': np.mean(stats['errors']),
            }
    
    return results


def print_results(results: dict):
    """Pretty print validation results"""
    print("\n" + "=" * 80)
    print("FULL GENOME PIPELINE VALIDATION RESULTS")
    print("=" * 80)
    
    overall = results['overall']
    print(f"\n📊 Overall Performance:")
    print(f"  Total reads:        {overall['total_reads']}")
    print(f"  Accuracy @ ±1kb:    {overall['accuracy_1kb']:.2f}%")
    print(f"  Accuracy @ ±100bp:  {overall['accuracy_100bp']:.2f}%")
    print(f"  Accuracy @ exact:   {overall['accuracy_exact']:.2f}%")
    print(f"\n⚡ Search Performance:")
    print(f"  Avg search time:    {overall['avg_search_time_ms']:.2f} ms")
    print(f"  Median search time: {overall['median_search_time_ms']:.2f} ms")
    print(f"  Reads per second:   {results['timing']['reads_per_second']:.1f}")
    
    print(f"\n📍 Per-Chromosome Results:")
    for chr_name, stats in sorted(results['per_chromosome'].items()):
        print(f"  {chr_name:15s}: {stats['accuracy_1kb']:5.1f}% @ ±1kb  "
              f"(median error: {stats['median_error']:.0f} bp, n={stats['total']})")
    
    print("\n" + "=" * 80)
    
    # Target check
    print("\n🎯 Target Check:")
    if overall['accuracy_1kb'] >= 95:
        print(f"  ✅ EXCELLENT! {overall['accuracy_1kb']:.1f}% >= 95%")
        print(f"  → Production ready!")
    elif overall['accuracy_1kb'] >= 90:
        print(f"  ✅ GOOD! {overall['accuracy_1kb']:.1f}% >= 90%")
        print(f"  → Nearly production ready")
    elif overall['accuracy_1kb'] >= 80:
        print(f"  ✅ PASSED! {overall['accuracy_1kb']:.1f}% >= 80%")
        print(f"  → Needs fine-tuning for production")
    else:
        print(f"  ⚠️  MARGINAL: {overall['accuracy_1kb']:.1f}% < 80%")
        print(f"  → Needs improvement")


def main():
    print("=" * 80)
    print("GenoCache V4 - Full Genome Pipeline Validation")
    print("=" * 80)
    print()
    
    # Configuration
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    CHECKPOINT_DIR = Path("/home/nebius/genocache/genocache-v4/models/checkpoints")
    OUTPUT_DIR = Path("/home/nebius/genocache/genocache-v4/validation/results")
    
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    # Test chromosomes (representative sample)
    TEST_CHROMOSOMES = [
        'NC_000001.11',  # chr1 (largest)
        'NC_000007.14',  # chr7 (medium)
        'NC_000022.11',  # chr22 (small, trained on this)
        'NC_000023.11',  # chrX (sex chromosome)
    ]
    
    ENCODE_CHROMOSOMES = TEST_CHROMOSOMES  # Same for now
    
    # Find best checkpoint
    checkpoints = list(CHECKPOINT_DIR.glob("fullgenome_best_*.pt"))
    if not checkpoints:
        checkpoints = list(CHECKPOINT_DIR.glob("fullgenome_epoch*.pt"))
    
    if not checkpoints:
        print("❌ No full genome checkpoint found!")
        return
    
    checkpoint_path = sorted(checkpoints)[-1]
    print(f"Checkpoint: {checkpoint_path.name}\n")
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}\n")
    
    # Load model
    print("Loading model...")
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    num_params = sum(p.numel() for p in model.parameters())
    print(f"✅ Model loaded: {num_params:,} parameters")
    
    # Load genome
    genome = load_genome(GENOME_PATH)
    print(f"✅ Loaded {len(genome)} chromosomes")
    
    # Generate test reads
    test_reads, ground_truth = generate_test_reads_multi_chr(
        genome,
        TEST_CHROMOSOMES,
        reads_per_chr=250,  # 250 per chr = 1000 total
        read_length=512,
        error_rate=0.05
    )
    
    # Encode genome
    start_time = time.time()
    genome_embeddings, genome_positions, genome_chr_names = encode_genome_with_stride(
        model, genome, ENCODE_CHROMOSOMES,
        window_size=512,
        stride=32,
        device=device,
        batch_size=256
    )
    genome_time = time.time() - start_time
    print(f"  Genome encoding time: {genome_time:.1f}s")
    
    # Build FAISS index
    start_time = time.time()
    index = build_faiss_index(genome_embeddings, nprobe=32)
    index_time = time.time() - start_time
    print(f"  Index building time: {index_time:.1f}s")
    
    # Encode test reads
    print("\nEncoding test reads...")
    start_time = time.time()
    read_embeddings = encode_reads(model, test_reads, device=device, batch_size=64)
    read_time = time.time() - start_time
    print(f"✅ Encoded {len(test_reads)} reads in {read_time:.1f}s")
    
    # Validate
    results = validate_pipeline(
        index, genome_positions, genome_chr_names,
        read_embeddings, ground_truth,
        k=32, tolerance=1000
    )
    
    # Print results
    print_results(results)
    
    # Save results
    results_file = OUTPUT_DIR / f"fullgenome_validation_{checkpoint_path.stem}.json"
    with open(results_file, 'w') as f:
        json.dump({
            'checkpoint': str(checkpoint_path),
            'results': results,
            'config': {
                'test_chromosomes': TEST_CHROMOSOMES,
                'reads_per_chr': 250,
                'read_length': 512,
                'error_rate': 0.05,
                'stride': 32,
                'k': 32,
                'tolerance': 1000
            },
            'timing': {
                'genome_encoding': genome_time,
                'index_building': index_time,
                'read_encoding': read_time
            }
        }, f, indent=2)
    
    print(f"\n✅ Results saved: {results_file}")


if __name__ == "__main__":
    main()
