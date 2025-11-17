#!/usr/bin/env python3
"""
Debug single read through entire pipeline to find position bug
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from seeding_nal import NALSeeding
from chaining_nal import NALChaining

# Read to debug: b095fff0-6e9d-4e30-92e8-d10ef3eac27e
# minimap2 maps to: Chr18:73838260
# NAL maps to: Chr2:239776577

def load_test_read():
    """Load the specific test read"""
    with open('data/test_100reads.fastq') as f:
        for line in f:
            if 'b095fff0-6e9d-4e30-92e8-d10ef3eac27e' in line:
                read_id = line.strip()[1:].split()[0]
                read_seq = f.readline().strip()
                return read_id, read_seq
    return None, None

def main():
    print("="*80)
    print("DEBUGGING SINGLE READ: b095fff0-6e9d-4e30-92e8-d10ef3eac27e")
    print("="*80)
    print()
    
    # Load read
    read_id, read_seq = load_test_read()
    if not read_seq:
        print("❌ Could not load read!")
        return
    
    print(f"Read ID: {read_id}")
    print(f"Read length: {len(read_seq)}")
    print(f"First 50bp: {read_seq[:50]}")
    print()
    
    print("Expected (minimap2):")
    print("  Chromosome: NC_000018.10 (Chr18)")
    print("  Position: 73,838,260")
    print()
    
    print("Actual (NAL):")
    print("  Chromosome: NC_000002.12 (Chr2)")
    print("  Position: 239,776,577")
    print()
    
    # Initialize seeder
    print("-"*80)
    print("Step 1: SEEDING")
    print("-"*80)
    
    seeder = NALSeeding(
        model_path='../models/genocache_nal.pt',
        index_path='../indexes/genocache_nal_stride32.index',
        positions_path='../indexes/genocache_nal_stride32.positions.npz',
        seed_len=512,
        K=32,
        device='cuda'
    )
    
    # Extract 6 seeds
    print(f"\nExtracting 6 seeds from read...")
    anchors = seeder.get_anchors(read_seq, num_seeds=6, K=32)
    
    print(f"\nGot {len(anchors)} anchors from FAISS")
    print("\nAnchor details:")
    
    # Group by seed
    from collections import defaultdict
    by_seed = defaultdict(list)
    for anchor in anchors:
        by_seed[anchor['seed_idx']].append(anchor)
    
    for seed_idx in sorted(by_seed.keys()):
        seed_anchors = by_seed[seed_idx]
        print(f"\n  Seed {seed_idx}: {len(seed_anchors)} anchors")
        
        # Count chromosomes
        chr_counts = defaultdict(int)
        for a in seed_anchors:
            chr_counts[a['ref_chr']] += 1
        
        print(f"    Chromosomes: {dict(chr_counts)}")
        
        # Show top 3 anchors by similarity
        top3 = sorted(seed_anchors, key=lambda x: x['similarity'], reverse=True)[:3]
        for i, a in enumerate(top3, 1):
            print(f"    #{i}: {a['ref_chr']}:{a['ref_pos']} (sim={a['similarity']:.3f})")
    
    # Check: Does Chr18 appear in anchors?
    chr18_anchors = [a for a in anchors if 'NC_000018' in a['ref_chr']]
    chr2_anchors = [a for a in anchors if 'NC_000002' in a['ref_chr']]
    
    print(f"\n📊 Chromosome distribution:")
    print(f"  Chr18 anchors: {len(chr18_anchors)}")
    print(f"  Chr2 anchors: {len(chr2_anchors)}")
    
    if chr18_anchors:
        print(f"\n  ✅ Chr18 IS in the anchor set!")
        print(f"     Positions: {[a['ref_pos'] for a in chr18_anchors[:5]]}")
    else:
        print(f"\n  ❌ Chr18 NOT in anchor set - problem in indexing/retrieval!")
    
    # Step 2: Chaining
    print("\n" + "-"*80)
    print("Step 2: CHAINING")
    print("-"*80)
    
    chainer = NALChaining(tolerance=1000, min_chain_score=3, top_k=5)
    chains, _ = chainer.chain_with_rescue_check(
        anchors, 
        read_len=len(read_seq), 
        num_seeds=6, 
        seed_len=512
    )
    
    print(f"\nGot {len(chains)} chains")
    
    for i, chain in enumerate(chains[:5], 1):
        print(f"\n  Chain {i}:")
        print(f"    Chr: {chain['ref_chr']}")
        print(f"    Pos: {chain['ref_pos']:,}")
        print(f"    Score: {chain['score']}")
        print(f"    Strand: {chain['strand']}")
        print(f"    Anchors: {len(chain.get('anchors', []))}")
    
    # Analysis
    print("\n" + "="*80)
    print("ANALYSIS")
    print("="*80)
    
    if not chr18_anchors:
        print("\n❌ ROOT CAUSE: Chr18 not in anchor set!")
        print("   Problem is in index building or FAISS retrieval")
        print("   The seeds from this read are not being matched to Chr18")
    elif len(chr18_anchors) < len(chr2_anchors):
        print(f"\n⚠️  Chr18 has fewer anchors ({len(chr18_anchors)}) than Chr2 ({len(chr2_anchors)})")
        print("   Chaining algorithm picks Chr2 because it has more anchors")
        print("   But minimap2 thinks Chr18 is correct!")
    else:
        print(f"\n❓ Chr18 has anchors but didn't win chaining")
        print("   Need to investigate chaining algorithm")
    
    # Check if any chain is on Chr18
    chr18_chains = [c for c in chains if 'NC_000018' in c['ref_chr']]
    if chr18_chains:
        print(f"\n   Chr18 chains exist (rank {chains.index(chr18_chains[0])+1})")
        print(f"   But Chr2 chain won with higher score")
    
    print()

if __name__ == '__main__':
    main()
