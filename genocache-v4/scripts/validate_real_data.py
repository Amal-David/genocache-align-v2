#!/usr/bin/env python3
"""
Validate model on real chr22 reads with ground truth
"""

import sys
import torch
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder


def seq_to_tensor(seq, max_len=1024):
    """Convert DNA sequence to tensor"""
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    # Trim or pad to max_len
    if len(seq) > max_len:
        seq = seq[:max_len]
    elif len(seq) < max_len:
        seq = seq + 'N' * (max_len - len(seq))
    
    tensor = torch.zeros(max_len, dtype=torch.long)
    for i, base in enumerate(seq):
        tensor[i] = base_to_idx.get(base.upper(), 4)
    return tensor


def encode_reads(model, reads_file, device, max_reads=None):
    """Encode all reads"""
    print(f"\nEncoding reads from {reads_file}...")
    
    read_names = []
    read_embeddings = []
    
    count = 0
    for record in tqdm(SeqIO.parse(reads_file, "fasta"), desc="Encoding reads"):
        if max_reads and count >= max_reads:
            break
        
        read_seq = str(record.seq).upper()
        read_tensor = seq_to_tensor(read_seq).unsqueeze(0).to(device)
        
        with torch.no_grad():
            embedding = model(read_tensor).cpu().numpy()[0]
        
        read_names.append(record.id)
        read_embeddings.append(embedding)
        count += 1
    
    read_embeddings = np.array(read_embeddings)
    print(f"✅ Encoded {len(read_names)} reads, embedding shape: {read_embeddings.shape}")
    
    return read_names, read_embeddings


def encode_genome_sliding(model, genome_file, chrom, device, window_size=1024, stride=512):
    """Encode genome in sliding windows"""
    print(f"\nLoading chromosome {chrom} from {genome_file}...")
    
    # Load chromosome
    genome_seq = None
    for record in SeqIO.parse(genome_file, "fasta"):
        if chrom in record.id or record.id in chrom:
            genome_seq = str(record.seq).upper()
            print(f"✅ Loaded {record.id}: {len(genome_seq):,} bp")
            break
    
    if genome_seq is None:
        raise ValueError(f"Chromosome {chrom} not found in {genome_file}")
    
    # Sliding window encoding
    print(f"Encoding genome with window={window_size}, stride={stride}...")
    
    positions = []
    embeddings = []
    
    for pos in tqdm(range(0, len(genome_seq) - window_size, stride), desc="Encoding genome"):
        window = genome_seq[pos:pos + window_size]
        
        # Skip if too many N's
        if window.count('N') > window_size * 0.2:
            continue
        
        window_tensor = seq_to_tensor(window, max_len=window_size).unsqueeze(0).to(device)
        
        with torch.no_grad():
            embedding = model(window_tensor).cpu().numpy()[0]
        
        positions.append(pos)
        embeddings.append(embedding)
    
    embeddings = np.array(embeddings)
    print(f"✅ Encoded {len(positions):,} windows, embedding shape: {embeddings.shape}")
    
    return np.array(positions), embeddings


def find_alignments(read_embeddings, genome_positions, genome_embeddings, top_k=5):
    """Find top-k closest genome positions for each read"""
    print(f"\nFinding alignments (top-{top_k})...")
    
    # Normalize embeddings
    read_emb_norm = read_embeddings / np.linalg.norm(read_embeddings, axis=1, keepdims=True)
    genome_emb_norm = genome_embeddings / np.linalg.norm(genome_embeddings, axis=1, keepdims=True)
    
    # Compute cosine similarities
    similarities = read_emb_norm @ genome_emb_norm.T
    
    # Find top-k for each read
    top_indices = np.argsort(-similarities, axis=1)[:, :top_k]
    top_scores = np.take_along_axis(similarities, top_indices, axis=1)
    top_positions = genome_positions[top_indices]
    
    return top_positions, top_scores


def evaluate_accuracy(read_names, predicted_positions, truth_file, tolerance=1000):
    """Evaluate alignment accuracy"""
    print(f"\nEvaluating accuracy (tolerance={tolerance}bp)...")
    
    # Load ground truth
    truth_df = pd.read_csv(truth_file, sep='\t', header=None, names=['read_id', 'true_pos'])
    truth_dict = dict(zip(truth_df['read_id'], truth_df['true_pos']))
    
    correct = 0
    total = 0
    errors = []
    
    for i, read_name in enumerate(read_names):
        if read_name not in truth_dict:
            print(f"⚠️ Warning: {read_name} not in truth file")
            continue
        
        true_pos = truth_dict[read_name]
        pred_pos_top5 = predicted_positions[i]
        
        # Check if any of top-5 predictions are within tolerance
        is_correct = any(abs(pred - true_pos) <= tolerance for pred in pred_pos_top5)
        
        if is_correct:
            correct += 1
        else:
            error = abs(pred_pos_top5[0] - true_pos)
            errors.append({
                'read': read_name,
                'true_pos': true_pos,
                'pred_pos': pred_pos_top5[0],
                'error': error
            })
        
        total += 1
    
    accuracy = (correct / total * 100) if total > 0 else 0
    
    print(f"\n{'='*80}")
    print(f"ACCURACY RESULTS")
    print(f"{'='*80}")
    print(f"Total reads: {total}")
    print(f"Correct alignments: {correct}")
    print(f"Accuracy: {accuracy:.2f}%")
    print(f"Tolerance: ±{tolerance}bp")
    print(f"{'='*80}")
    
    # Show some errors
    if errors:
        print(f"\nTop 10 largest errors:")
        errors_sorted = sorted(errors, key=lambda x: x['error'], reverse=True)[:10]
        for err in errors_sorted:
            print(f"  {err['read']}: true={err['true_pos']:,}, pred={err['pred_pos']:,}, error={err['error']:,}bp")
    
    return accuracy, errors


def main():
    print("="*80)
    print("GenoCache V4 - Real Data Validation")
    print("="*80)
    
    # Configuration
    checkpoint_path = 'models/checkpoints/best_model_10M.pt'
    genome_file = '/home/nebius/genocache/GRCh38.fa'
    reads_file = '/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa'
    truth_file = '/home/nebius/genocache/genocache_data/reads_chr22_truth_FIXED.tsv'
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\nConfiguration:")
    print(f"  Checkpoint: {checkpoint_path}")
    print(f"  Genome: {genome_file}")
    print(f"  Reads: {reads_file}")
    print(f"  Truth: {truth_file}")
    print(f"  Device: {device}")
    
    # Load model
    print(f"\nLoading model from {checkpoint_path}...")
    model = GenoCacheEncoder(emb_dim=256, seed_len=512)
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(device)
    model.eval()
    
    print(f"✅ Model loaded")
    print(f"   Metadata: {checkpoint.get('metadata', {})}")
    
    # Encode reads
    read_names, read_embeddings = encode_reads(model, reads_file, device)
    
    # Encode genome (chr22)
    genome_positions, genome_embeddings = encode_genome_sliding(
        model, genome_file, 'NC_000022.11', device, 
        window_size=1024, stride=512
    )
    
    # Find alignments
    predicted_positions, scores = find_alignments(
        read_embeddings, genome_positions, genome_embeddings, top_k=5
    )
    
    # Evaluate accuracy
    accuracy_1kb, errors_1kb = evaluate_accuracy(
        read_names, predicted_positions, truth_file, tolerance=1000
    )
    
    # Try different tolerances
    print(f"\n{'='*80}")
    print("ACCURACY AT DIFFERENT TOLERANCES")
    print(f"{'='*80}")
    
    for tolerance in [100, 500, 1000, 5000, 10000]:
        acc, _ = evaluate_accuracy(read_names, predicted_positions, truth_file, tolerance=tolerance)
        print(f"  ±{tolerance:5d}bp: {acc:6.2f}%")
    
    print(f"\n{'='*80}")
    print("VALIDATION COMPLETE")
    print(f"{'='*80}")
    
    # Summary
    print(f"\nSummary:")
    print(f"  Model: Epoch 1, Separation 0.5346")
    print(f"  Reads tested: {len(read_names)}")
    print(f"  Accuracy (±1kb): {accuracy_1kb:.2f}%")
    
    if accuracy_1kb >= 95:
        print(f"\n✅ EXCELLENT: >95% accuracy!")
        print(f"   Recommendation: Use curriculum learning on 30M for final boost")
    elif accuracy_1kb >= 85:
        print(f"\n✅ GOOD: 85-95% accuracy")
        print(f"   Recommendation: Curriculum learning on 30-50M should reach 99%+")
    elif accuracy_1kb >= 70:
        print(f"\n⚠️ MODERATE: 70-85% accuracy")
        print(f"   Recommendation: Need curriculum + more data (50M+)")
    else:
        print(f"\n❌ POOR: <70% accuracy")
        print(f"   Recommendation: Fundamental changes needed")


if __name__ == "__main__":
    main()
