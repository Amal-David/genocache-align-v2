"""
Validation script for GenoCache-Align V2
Tests model recall on synthetic queries
"""

import torch
import faiss
import numpy as np
from pathlib import Path
import json
from tqdm import tqdm
import argparse

from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode


def generate_test_queries(fasta_path, chrom, num_queries=1000, read_len=1000):
    """Generate random test queries from reference"""
    print(f"Generating {num_queries} test queries...")
    
    ref = ReferenceGenome(fasta_path)
    seq = ref.sequences[chrom]
    seq_len = len(seq)
    
    queries = []
    positions = []
    
    for _ in range(num_queries):
        # Random position
        pos = np.random.randint(0, seq_len - read_len)
        query_seq = seq[pos:pos + read_len]
        
        # Skip if contains too many Ns
        if query_seq.count('N') > read_len * 0.1:
            continue
        
        queries.append(query_seq)
        positions.append(pos)
    
    print(f"✅ Generated {len(queries)} queries")
    return queries, positions


def validate_model(checkpoint_path, index_path, fasta_path, chrom, 
                  num_queries=1000, seed_offsets=[0, 128, 256, 384, 512]):
    """
    Validate model accuracy
    
    Args:
        checkpoint_path: Path to trained model
        index_path: Path to FAISS index
        fasta_path: Reference FASTA
        chrom: Chromosome
        num_queries: Number of test queries
        seed_offsets: Positions to extract seeds from reads
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n{'='*80}")
    print("VALIDATION")
    print(f"{'='*80}")
    print(f"Checkpoint: {checkpoint_path}")
    print(f"Index: {index_path}")
    print(f"Device: {device}")
    print(f"{'='*80}\n")
    
    # Load model
    print("Loading model...")
    model, metadata = ImprovedCNN.load_checkpoint(checkpoint_path, device=device)
    model.eval()
    
    seed_len = model.input_len
    
    # Load index
    print("Loading index...")
    index = faiss.read_index(str(index_path))
    print(f"✅ Index: {index.ntotal:,} vectors")
    
    # Load positions
    manifest_dir = Path(index_path).parent
    manifest_path = manifest_dir / 'indexes_manifest.json'
    
    with open(manifest_path) as f:
        manifest = json.load(f)
    
    positions_path = manifest['source']['reference']['positions']['path']
    ref_positions = np.load(positions_path)
    ref_stride = manifest['source']['reference']['stride']
    
    print(f"✅ Positions: {len(ref_positions):,}")
    print(f"   Range: {ref_positions.min():,} - {ref_positions.max():,}")
    print(f"   Stride: {ref_stride}")
    
    # Generate test queries
    queries, true_positions = generate_test_queries(
        fasta_path, chrom, num_queries, read_len=1024
    )
    
    # Validate
    print(f"\nValidating with {len(seed_offsets)} seeds per read...")
    
    results = {
        'top1': 0,
        'top5': 0,
        'top10': 0,
        'top32': 0,
        'total': 0
    }
    
    with torch.no_grad():
        for query_seq, true_pos in tqdm(zip(queries, true_positions), total=len(queries)):
            # Extract seeds at multiple offsets
            query_found = False
            
            for offset in seed_offsets:
                if offset + seed_len > len(query_seq):
                    continue
                
                seed_seq = query_seq[offset:offset + seed_len]
                seed_tensor = one_hot_encode(seed_seq).unsqueeze(0).to(device)
                
                # Encode
                embedding = model(seed_tensor).cpu().numpy()
                
                # Search
                D, I = index.search(embedding, k=32)
                
                # Get predicted positions
                pred_positions = ref_positions[I[0]]
                
                # True seed position in reference
                true_seed_pos = true_pos + offset
                
                # Check if any prediction is within tolerance
                tolerance = ref_stride * 2  # Allow ±2 strides
                distances = np.abs(pred_positions - true_seed_pos)
                
                if distances[0] <= tolerance:
                    results['top1'] += 1
                    query_found = True
                    break
                elif distances[:5].min() <= tolerance:
                    results['top5'] += 1
                    query_found = True
                    break
                elif distances[:10].min() <= tolerance:
                    results['top10'] += 1
                    query_found = True
                    break
                elif distances[:32].min() <= tolerance:
                    results['top32'] += 1
                    query_found = True
                    break
            
            if query_found:
                results['total'] += 1
    
    # Calculate accuracies
    num_queries = len(queries)
    
    print(f"\n{'='*80}")
    print("RESULTS")
    print(f"{'='*80}")
    print(f"Test queries: {num_queries}")
    print(f"Seeds per read: {len(seed_offsets)}")
    print(f"Tolerance: ±{ref_stride * 2} bp")
    print()
    print(f"Top-1 accuracy:  {results['top1']:4d} / {num_queries} = {results['top1']/num_queries*100:.2f}%")
    print(f"Top-5 accuracy:  {results['top5']:4d} / {num_queries} = {results['top5']/num_queries*100:.2f}%")
    print(f"Top-10 accuracy: {results['top10']:4d} / {num_queries} = {results['top10']/num_queries*100:.2f}%")
    print(f"Top-32 accuracy: {results['top32']:4d} / {num_queries} = {results['top32']/num_queries*100:.2f}%")
    print(f"Total recall:    {results['total']:4d} / {num_queries} = {results['total']/num_queries*100:.2f}%")
    print(f"{'='*80}\n")
    
    # Expected performance
    print("Expected performance (based on h100 validation):")
    print("  Top-1:  >60%")
    print("  Top-10: >85%")
    print("  Top-32: >90%")
    
    if results['top1'] / num_queries < 0.5:
        print("\n⚠️  WARNING: Top-1 accuracy below 50%")
        print("   Possible issues:")
        print("   - Model underfitted (train longer)")
        print("   - Embeddings not normalized")
        print("   - Wrong checkpoint used for encoding")
    elif results['top1'] / num_queries > 0.6:
        print("\n✅ Performance looks good!")
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate GenoCache-Align V2')
    
    parser.add_argument('--checkpoint', type=str, required=True,
                       help='Path to model checkpoint')
    parser.add_argument('--index', type=str, required=True,
                       help='Path to FAISS index')
    parser.add_argument('--fasta', type=str, required=True,
                       help='Path to reference FASTA')
    parser.add_argument('--chrom', type=str, required=True,
                       help='Chromosome to test on')
    parser.add_argument('--num-queries', type=int, default=1000,
                       help='Number of test queries')
    
    args = parser.parse_args()
    
    validate_model(
        checkpoint_path=args.checkpoint,
        index_path=args.index,
        fasta_path=args.fasta,
        chrom=args.chrom,
        num_queries=args.num_queries
    )
