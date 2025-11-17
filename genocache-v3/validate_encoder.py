#!/usr/bin/env python3
"""
Validation script for GenoCacheEncoder
Tests inference speed and embedding quality
"""

import os
import sys
import time
import random
import argparse
import numpy as np
import torch
import torch.nn.functional as F
from Bio import SeqIO

from genocache_encoder import GenoCacheEncoder, seq_to_tokens, augment_sequence


def test_inference_speed(model, device, batch_size=256, n_batches=100):
    """Test inference speed (throughput and latency)"""
    print("="*80)
    print("INFERENCE SPEED TEST")
    print("="*80)
    
    model.eval()
    
    # Generate random test sequences
    test_sequences = torch.randint(0, 4, (batch_size, 512), device=device)
    
    # Warmup
    with torch.no_grad():
        for _ in range(10):
            _ = model(test_sequences)
    
    # Measure throughput
    torch.cuda.synchronize() if device == "cuda" else None
    start = time.time()
    
    with torch.no_grad():
        for _ in range(n_batches):
            _ = model(test_sequences)
    
    torch.cuda.synchronize() if device == "cuda" else None
    end = time.time()
    
    total_seqs = batch_size * n_batches
    total_time = end - start
    throughput = total_seqs / total_time
    latency_per_seed = (total_time / total_seqs) * 1e6  # microseconds
    
    print(f"Batch size: {batch_size}")
    print(f"Total sequences: {total_seqs:,}")
    print(f"Total time: {total_time:.2f}s")
    print(f"Throughput: {throughput:,.0f} seeds/sec")
    print(f"Latency per seed: {latency_per_seed:.1f} μs")
    print()
    
    # Target check
    target_latency = 25  # μs
    if latency_per_seed < target_latency:
        print(f"✅ PASS: {latency_per_seed:.1f} μs < {target_latency} μs target")
    else:
        print(f"⚠️  MISS: {latency_per_seed:.1f} μs > {target_latency} μs target")
    
    print()
    return throughput, latency_per_seed


def test_embedding_quality(model, device, ref_sequence, n_tests=100):
    """
    Test embedding quality:
    - Similar sequences should have high cosine similarity
    - Different sequences should have low cosine similarity
    """
    print("="*80)
    print("EMBEDDING QUALITY TEST")
    print("="*80)
    
    model.eval()
    
    seed_len = 512
    similarities_same = []
    similarities_diff = []
    
    with torch.no_grad():
        for _ in range(n_tests):
            # Sample a random sequence
            pos = random.randint(0, len(ref_sequence) - seed_len)
            seq = ref_sequence[pos:pos + seed_len]
            
            # Create augmented version (should be similar)
            aug_seq = augment_sequence(seq, error_rate_range=(0.02, 0.05))
            
            # Sample a different sequence (should be different)
            diff_pos = random.randint(0, len(ref_sequence) - seed_len)
            # Ensure it's far away
            while abs(diff_pos - pos) < 10000:
                diff_pos = random.randint(0, len(ref_sequence) - seed_len)
            diff_seq = ref_sequence[diff_pos:diff_pos + seed_len]
            
            # Encode
            seq_tokens = torch.tensor([seq_to_tokens(seq)], device=device)
            aug_tokens = torch.tensor([seq_to_tokens(aug_seq)], device=device)
            diff_tokens = torch.tensor([seq_to_tokens(diff_seq)], device=device)
            
            seq_emb = model(seq_tokens)
            aug_emb = model(aug_tokens)
            diff_emb = model(diff_tokens)
            
            # Compute similarities
            sim_same = (seq_emb * aug_emb).sum().item()
            sim_diff = (seq_emb * diff_emb).sum().item()
            
            similarities_same.append(sim_same)
            similarities_diff.append(sim_diff)
    
    # Statistics
    mean_same = np.mean(similarities_same)
    std_same = np.std(similarities_same)
    mean_diff = np.mean(similarities_diff)
    std_diff = np.std(similarities_diff)
    
    print(f"Same sequence (augmented):")
    print(f"  Mean similarity: {mean_same:.3f} ± {std_same:.3f}")
    print()
    print(f"Different sequences:")
    print(f"  Mean similarity: {mean_diff:.3f} ± {std_diff:.3f}")
    print()
    
    # Target check
    target_same = 0.7
    target_diff = 0.3
    
    if mean_same > target_same and mean_diff < target_diff:
        print(f"✅ PASS: Good separation ({mean_same:.3f} > {target_same}, {mean_diff:.3f} < {target_diff})")
    else:
        print(f"⚠️  MISS: Insufficient separation")
    
    print()
    return mean_same, mean_diff


def test_translation_continuity(model, device, ref_sequence, n_tests=50):
    """
    Test translation continuity:
    - Embeddings should change smoothly for shifted sequences
    """
    print("="*80)
    print("TRANSLATION CONTINUITY TEST")
    print("="*80)
    
    model.eval()
    
    seed_len = 512
    shifts = list(range(0, 257, 8))  # Test shifts from 0 to 256
    
    avg_similarities = []
    
    with torch.no_grad():
        for _ in range(n_tests):
            # Sample a random sequence
            pos = random.randint(0, len(ref_sequence) - seed_len - 256)
            long_seq = ref_sequence[pos:pos + seed_len + 256]
            
            # Get embeddings for different shifts
            base_seq = long_seq[:seed_len]
            base_tokens = torch.tensor([seq_to_tokens(base_seq)], device=device)
            base_emb = model(base_tokens)
            
            similarities = []
            for shift in shifts[1:]:  # Skip 0 (would be 1.0)
                shifted_seq = long_seq[shift:shift + seed_len]
                shifted_tokens = torch.tensor([seq_to_tokens(shifted_seq)], device=device)
                shifted_emb = model(shifted_tokens)
                
                sim = (base_emb * shifted_emb).sum().item()
                similarities.append(sim)
            
            avg_similarities.append(similarities)
    
    # Average across tests
    avg_similarities = np.mean(avg_similarities, axis=0)
    
    print(f"Shift (bp) | Similarity")
    print("-" * 30)
    for shift, sim in zip(shifts[1:], avg_similarities):
        print(f"{shift:>8}   | {sim:.3f}")
    
    print()
    
    # Check that similarity decreases smoothly
    # At shift=seed_len/2 (256), similarity should be moderate
    mid_idx = len(shifts) // 2
    mid_sim = avg_similarities[mid_idx]
    
    if mid_sim > 0.3 and mid_sim < 0.7:
        print(f"✅ PASS: Good translation continuity (mid-shift sim = {mid_sim:.3f})")
    else:
        print(f"⚠️  CHECK: Translation continuity (mid-shift sim = {mid_sim:.3f})")
    
    print()


def main():
    parser = argparse.ArgumentParser(description="Validate GenoCacheEncoder")
    parser.add_argument("--model", type=str, required=True,
                       help="Path to model checkpoint")
    parser.add_argument("--fasta", type=str, default="../GRCh38.fa",
                       help="Reference FASTA")
    parser.add_argument("--chrom-id", type=str, default="NC_000001.11",
                       help="Chromosome ID for testing")
    parser.add_argument("--seed-len", type=int, default=512,
                       help="Seed length")
    parser.add_argument("--emb-dim", type=int, default=256,
                       help="Embedding dimension")
    
    args = parser.parse_args()
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    print("="*80)
    print("GenoCache V3 Validation")
    print("="*80)
    print(f"Model: {args.model}")
    print(f"Device: {device}")
    print()
    
    # Load model
    print("Loading model...")
    model = GenoCacheEncoder(seed_len=args.seed_len, emb_dim=args.emb_dim)
    
    checkpoint = torch.load(args.model, map_location=device)
    if 'model_state_dict' in checkpoint:
        model.load_state_dict(checkpoint['model_state_dict'])
        print(f"  Epoch: {checkpoint.get('epoch', 'unknown')}")
        print(f"  Loss: {checkpoint.get('loss', 'unknown'):.6f}")
    else:
        model.load_state_dict(checkpoint)
    
    model = model.to(device)
    model.eval()
    print("✅ Model loaded")
    print()
    
    # Load reference sequence
    print(f"Loading reference: {args.chrom_id}...")
    ref_sequence = None
    for record in SeqIO.parse(args.fasta, "fasta"):
        if record.id == args.chrom_id:
            ref_sequence = str(record.seq).upper()
            break
    
    if ref_sequence is None:
        print(f"❌ Could not find {args.chrom_id} in {args.fasta}")
        return
    
    print(f"✅ Loaded {len(ref_sequence):,} bp")
    print()
    
    # Run tests
    test_inference_speed(model, device)
    test_embedding_quality(model, device, ref_sequence)
    test_translation_continuity(model, device, ref_sequence)
    
    print("="*80)
    print("VALIDATION COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
