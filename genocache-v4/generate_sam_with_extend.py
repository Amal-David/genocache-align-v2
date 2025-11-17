#!/usr/bin/env python3
"""
Generate SAM Output with EXTEND Phase

This produces minimap2-style SAM output with:
- CIGAR strings
- MAPQ scores
- Proper SAM format
- Using EXTEND phase (THE FIX!)

Expected: 37% → 95%+ chromosome accuracy compared to OLD method
"""

import sys
import time
from pathlib import Path
sys.path.append(str(Path(__file__).parent))
sys.path.append(str(Path(__file__).parent / "models"))

import torch
import pickle
import faiss
from models.encoder import GenoCacheEncoder
from extend_phase import ExtendPhase
from fast_alignment import FastAligner
import parasail

def load_fasta(fasta_file):
    """Load FASTA file"""
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences.append((current_id, ''.join(current_seq)))
    
    return sequences

def write_sam_header(sam_file):
    """Write SAM header"""
    with open(sam_file, 'w') as f:
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        # Add main chromosomes (simplified)
        chrs = [
            ("NC_000022.11", 50818468),
            ("NC_000013.11", 114364328),
            ("NC_000014.9", 107043718),
            ("NC_000016.10", 90338345),
        ]
        for chr_name, length in chrs:
            f.write(f"@SQ\tSN:{chr_name}\tLN:{length}\n")
        f.write("@PG\tID:genocache\tPN:genocache-extend\tVN:4.0\tCL:with_extend_phase\n")

def calculate_mapq(score, second_best=0):
    """Calculate MAPQ from alignment score"""
    if second_best == 0 or score <= second_best:
        return 60
    score_diff = score - second_best
    if score_diff > 100:
        return 60
    elif score_diff > 50:
        return 40
    elif score_diff > 20:
        return 20
    else:
        return 10

def write_sam_alignment(sam_file, read_id, read_seq, alignment):
    """Write alignment to SAM file"""
    with open(sam_file, 'a') as f:
        if alignment is None:
            # Unmapped
            f.write(f"{read_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
        else:
            # Mapped
            flag = 0
            chr_name = alignment['chr']
            pos = alignment['start'] + 1  # SAM is 1-based
            mapq = alignment.get('mapq', 60)
            cigar = alignment.get('cigar', f"{len(read_seq)}M")
            
            # Tags
            tags = []
            tags.append(f"AS:i:{alignment.get('alignment_score', 0)}")
            tags.append(f"XS:i:{alignment.get('num_seeds', 0)}")
            tags.append(f"XC:i:{alignment.get('candidates_tested', 0)}")
            tags_str = '\t'.join(tags)
            
            f.write(f"{read_id}\t{flag}\t{chr_name}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{read_seq}\t*\t{tags_str}\n")

def simple_seeding_top_k(model, index, metadata, read_seq, top_k=5):
    """Simple seeding that returns top-k candidates"""
    # Encode read
    device = next(model.parameters()).device
    
    # Simple encoding (this is a placeholder - real implementation would chunk the read)
    # For now, just encode first 512bp
    seq_chunk = read_seq[:512].upper()
    char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    indices = [char_to_idx.get(c, 4) for c in seq_chunk]
    
    # Pad to 512
    while len(indices) < 512:
        indices.append(4)
    indices = indices[:512]
    
    # Encode
    x = torch.tensor([indices], dtype=torch.long).to(device)
    with torch.no_grad():
        embedding = model(x)
        embedding = embedding.cpu().numpy()[0]
    
    # Search FAISS
    D, I = index.search(embedding.reshape(1, -1).astype('float32'), top_k * 10)
    
    # Group by chromosome and get top-k chromosomes
    chr_candidates = {}
    for dist, idx in zip(D[0], I[0]):
        if idx < 0 or idx >= len(metadata['positions']):
            continue
        
        pos = metadata['positions'][idx]
        chr_name = metadata['chr_names'][idx]
        
        if chr_name not in chr_candidates:
            chr_candidates[chr_name] = {
                'chr': chr_name,
                'positions': [],
                'scores': [],
                'num_seeds': 0
            }
        
        chr_candidates[chr_name]['positions'].append(pos)
        chr_candidates[chr_name]['scores'].append(float(1.0 / (1.0 + dist)))
        chr_candidates[chr_name]['num_seeds'] += 1
    
    # Get top-k by seed count
    sorted_chrs = sorted(chr_candidates.values(), key=lambda x: x['num_seeds'], reverse=True)[:top_k]
    
    # Create candidate regions
    candidates = []
    for cand in sorted_chrs:
        if not cand['positions']:
            continue
        
        start = min(cand['positions'])
        end = max(cand['positions']) + 512
        
        candidates.append({
            'chr': cand['chr'],
            'start': start,
            'end': end,
            'num_seeds': cand['num_seeds'],
            'score': sum(cand['scores']) / len(cand['scores']) if cand['scores'] else 0
        })
    
    return candidates if candidates else None

def main():
    print("=" * 80)
    print("GENOCACHE WITH EXTEND PHASE - SAM OUTPUT GENERATION")
    print("=" * 80)
    print()
    
    # Configuration
    model_path = 'models/checkpoints/fullgenome_best_sep12.4035_epoch30.pt'
    index_path = 'indexes/genocache_v4_production.index'
    metadata_path = 'indexes/genocache_v4_production.metadata.pkl'
    test_reads = 'test_10_reads_exact.fa'
    output_sam = 'genocache_with_extend.sam'
    ref_genome = '../GRCh38.fa'
    
    # Check files
    for path, name in [(model_path, 'Model'), (index_path, 'Index'), 
                       (metadata_path, 'Metadata'), (test_reads, 'Test reads'),
                       (ref_genome, 'Reference')]:
        if not Path(path).exists():
            print(f"❌ {name} not found: {path}")
            return 1
    
    print("Loading components...")
    print()
    
    # Load model
    print("1. Loading model...")
    checkpoint = torch.load(model_path, map_location='cuda' if torch.cuda.is_available() else 'cpu', weights_only=True)
    model = GenoCacheEncoder(emb_dim=128, seed_len=512)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.cuda() if torch.cuda.is_available() else model
    model.eval()
    print(f"   ✅ Model loaded")
    
    # Load index
    print("2. Loading FAISS index...")
    index = faiss.read_index(index_path)
    print(f"   ✅ Index loaded ({index.ntotal:,} vectors)")
    
    # Load metadata
    print("3. Loading metadata...")
    with open(metadata_path, 'rb') as f:
        metadata = pickle.load(f)
    print(f"   ✅ Metadata loaded")
    
    # Load reference genome into memory (main chromosomes only)
    print("4. Loading reference genome...")
    genome_dict = {}
    current_chr = None
    current_seq = []
    target_chrs = ['NC_000022.11', 'NC_000013.11', 'NC_000014.9', 'NC_000016.10',
                   'NC_000019.10', 'NC_000001.11', 'NC_000021.9']  # Most common
    
    with open(ref_genome) as f:
        for line in f:
            if line.startswith('>'):
                if current_chr and current_chr in target_chrs:
                    genome_dict[current_chr] = ''.join(current_seq)
                    print(f"   Loaded: {current_chr} ({len(genome_dict[current_chr]):,} bp)")
                
                current_chr = line[1:].split()[0]
                current_seq = []
            elif current_chr in target_chrs:
                current_seq.append(line.strip())
        
        if current_chr and current_chr in target_chrs:
            genome_dict[current_chr] = ''.join(current_seq)
            print(f"   Loaded: {current_chr} ({len(genome_dict[current_chr]):,} bp)")
    
    print(f"   ✅ Loaded {len(genome_dict)} chromosomes")
    
    # Initialize aligner
    print("5. Initializing fast aligner...")
    aligner = FastAligner(genome_dict, mode='semi-global')
    print(f"   ✅ Aligner ready (parasail)")
    
    # Initialize EXTEND phase
    print("6. Initializing EXTEND phase...")
    extend = ExtendPhase(aligner, min_score_threshold=100)
    print(f"   ✅ EXTEND phase ready")
    print()
    
    # Load reads
    print("Loading test reads...")
    reads = load_fasta(test_reads)
    print(f"  ✅ Loaded {len(reads)} reads")
    print()
    
    # Initialize SAM
    print(f"Initializing SAM output: {output_sam}")
    write_sam_header(output_sam)
    print(f"  ✅ SAM header written")
    print()
    
    print("=" * 80)
    print("RUNNING PIPELINE WITH EXTEND PHASE")
    print("=" * 80)
    print()
    
    total_time = 0
    results = []
    
    for i, (read_id, read_seq) in enumerate(reads):
        print(f"Read {i+1}/{len(reads)}: {read_id} ({len(read_seq)} bp)")
        
        start_time = time.time()
        
        try:
            # Step 1: Seeding (get top-5 candidates)
            print(f"  [1/3] Seeding...", end=" ", flush=True)
            candidates = simple_seeding_top_k(model, index, metadata, read_seq, top_k=5)
            
            if not candidates:
                print("no candidates")
                write_sam_alignment(output_sam, read_id, read_seq, None)
                print(f"  ❌ Unmapped (no seeds)")
                continue
            
            print(f"found {len(candidates)} candidates")
            for j, cand in enumerate(candidates[:3]):
                print(f"     {j+1}. {cand['chr']} ({cand['num_seeds']} seeds)")
            
            # Step 2: EXTEND phase
            print(f"  [2/3] EXTEND (align to each)...", end=" ", flush=True)
            alignment = extend.extend_and_score(read_seq, candidates)
            
            if not alignment:
                print("no valid alignment")
                write_sam_alignment(output_sam, read_id, read_seq, None)
                print(f"  ❌ Unmapped (low score)")
                continue
            
            print(f"best: {alignment['chr']}")
            print(f"     Position: {alignment['start']}:{alignment['end']}")
            print(f"     Score: {alignment['alignment_score']}")
            print(f"     Status: {alignment['status']}")
            
            # Step 3: Write to SAM
            print(f"  [3/3] Writing to SAM...", end=" ", flush=True)
            write_sam_alignment(output_sam, read_id, read_seq, alignment)
            print("done")
            
            elapsed = time.time() - start_time
            total_time += elapsed
            print(f"  Time: {elapsed:.3f}s")
            print()
            
            results.append({'read_id': read_id, 'mapped': True})
            
        except Exception as e:
            print(f"\n  ❌ Error: {e}")
            import traceback
            traceback.print_exc()
            write_sam_alignment(output_sam, read_id, read_seq, None)
            print()
    
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    
    mapped = sum(1 for r in results if r.get('mapped'))
    print(f"Total reads: {len(results)}")
    print(f"Mapped: {mapped}/{len(results)} ({100*mapped/len(results):.1f}%)")
    print(f"Total time: {total_time:.3f}s")
    print(f"Avg time: {total_time/len(results):.3f}s per read")
    print()
    
    print(f"✅ SAM output: {output_sam}")
    print()
    print("Next steps:")
    print(f"1. Compare with OLD method:")
    print(f"   python3 compare_chromosome_accuracy.py {output_sam} minimap2_same_10reads.sam")
    print()
    print(f"2. View SAM file:")
    print(f"   head -30 {output_sam}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
