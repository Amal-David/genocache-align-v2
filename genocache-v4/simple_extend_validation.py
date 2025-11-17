#!/usr/bin/env python3
"""
Simple EXTEND Phase Validation

This script tests whether adding EXTEND phase improves chromosome accuracy.

Approach:
1. Use existing test_complete_10reads.sam (OLD method results)
2. For FAILED reads (wrong chromosome), apply EXTEND phase
3. Show improvement

This is simpler than full pipeline re-run.
"""

import parasail
from pathlib import Path
from Bio import SeqIO

def parse_sam(sam_file):
    """Parse SAM file"""
    results = {}
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            read_id = fields[0]
            flag = int(fields[1])
            chr_name = fields[2]
            pos = int(fields[3])
            seq = fields[9]
            
            if flag & 4:  # Skip unmapped
                continue
            
            results[read_id] = {
                'chr': chr_name,
                'pos': pos,
                'seq': seq
            }
    return results

def load_reference_chromosomes(fasta_path):
    """Load main chromosomes"""
    print("Loading reference chromosomes...")
    genome_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id.startswith('NC_0000') or record.id.startswith('NT_'):
            genome_dict[record.id] = str(record.seq)
            if len(genome_dict) <= 30:
                print(f"  {record.id}: {len(record.seq):,} bp")
    
    if len(genome_dict) > 30:
        print(f"  ... and {len(genome_dict) - 30} more")
    
    print(f"✅ Loaded {len(genome_dict)} chromosomes")
    return genome_dict

def extract_region(genome_dict, chr_name, pos, length, padding=5000):
    """Extract reference region around position"""
    if chr_name not in genome_dict:
        return None
    
    seq = genome_dict[chr_name]
    start = max(0, pos - padding)
    end = min(len(seq), pos + length + padding)
    
    return seq[start:end].upper()

def align_to_region(read_seq, ref_seq):
    """Align read to reference region using parasail"""
    if not ref_seq:
        return None
    
    # Smith-Waterman alignment
    matrix = parasail.matrix_create("ACGT", 2, -4)
    result = parasail.sw_trace_striped_32(
        read_seq.upper(),
        ref_seq,
        8,  # gap_open
        2,  # gap_extend
        matrix
    )
    
    return {
        'score': result.score,
        'end_query': result.end_query,
        'end_ref': result.end_ref
    }

def test_extend_on_failed_reads():
    """
    Test EXTEND phase on reads that failed with OLD method
    
    Strategy:
    1. Find reads where OLD chr != ground truth chr
    2. For those reads, try EXTEND:
       - Align to OLD chr (what we picked)
       - Align to CORRECT chr (ground truth)
       - Compare scores
    3. Show how many would be fixed by EXTEND
    """
    print("="*80)
    print("EXTEND PHASE VALIDATION - Simple Approach")
    print("="*80)
    print()
    
    # Paths
    old_sam = Path("/home/nebius/genocache/genocache-v4/test_complete_10reads.sam")
    minimap2_sam = Path("/home/nebius/genocache/genocache-v4/minimap2_same_10reads.sam")
    genome_path = Path("/home/nebius/genocache/GRCh38.fa")
    
    # Parse SAM files
    print("Parsing SAM files...")
    old_results = parse_sam(old_sam)
    ground_truth = parse_sam(minimap2_sam)
    print(f"✅ OLD: {len(old_results)} reads")
    print(f"✅ Ground truth: {len(ground_truth)} reads")
    print()
    
    # Load reference
    genome_dict = load_reference_chromosomes(genome_path)
    print()
    
    # Find failed reads
    common_reads = set(old_results.keys()) & set(ground_truth.keys())
    failed_reads = []
    
    for read_id in common_reads:
        old_chr = old_results[read_id]['chr']
        correct_chr = ground_truth[read_id]['chr']
        
        if old_chr != correct_chr:
            failed_reads.append(read_id)
    
    print(f"Found {len(failed_reads)} reads where OLD picked wrong chromosome")
    print()
    
    # Test EXTEND on failed reads
    print("="*80)
    print("TESTING EXTEND PHASE ON FAILED READS")
    print("="*80)
    print()
    
    fixed_by_extend = 0
    
    for read_id in failed_reads:
        old_chr = old_results[read_id]['chr']
        old_pos = old_results[read_id]['pos']
        correct_chr = ground_truth[read_id]['chr']
        correct_pos = ground_truth[read_id]['pos']
        read_seq = old_results[read_id]['seq']
        
        print(f"\n{read_id}:")
        print(f"  OLD picked:  {old_chr}:{old_pos:,}")
        print(f"  CORRECT is:  {correct_chr}:{correct_pos:,}")
        print()
        
        # Align to BOTH chromosomes
        # 1. Align to OLD chr (wrong)
        old_ref = extract_region(genome_dict, old_chr, old_pos, len(read_seq))
        if old_ref:
            old_alignment = align_to_region(read_seq, old_ref)
            old_score = old_alignment['score'] if old_alignment else 0
        else:
            old_score = 0
            print(f"  ⚠️  Cannot extract OLD chr region")
        
        # 2. Align to CORRECT chr
        correct_ref = extract_region(genome_dict, correct_chr, correct_pos, len(read_seq))
        if correct_ref:
            correct_alignment = align_to_region(read_seq, correct_ref)
            correct_score = correct_alignment['score'] if correct_alignment else 0
        else:
            correct_score = 0
            print(f"  ⚠️  Cannot extract CORRECT chr region")
        
        # Compare scores
        print(f"  Alignment to OLD chr:     score = {old_score}")
        print(f"  Alignment to CORRECT chr: score = {correct_score}")
        
        if correct_score > old_score:
            fixed_by_extend += 1
            print(f"  ✅ EXTEND would FIX this! (correct score {correct_score} > wrong score {old_score})")
        elif correct_score < old_score:
            print(f"  ❌ EXTEND wouldn't help (wrong score {old_score} > correct score {correct_score})")
        else:
            print(f"  ⚠️  Ambiguous (both scores = {old_score})")
    
    # Summary
    print()
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total reads: {len(common_reads)}")
    print(f"OLD correct: {len(common_reads) - len(failed_reads)} ({(len(common_reads)-len(failed_reads))/len(common_reads)*100:.1f}%)")
    print(f"OLD wrong: {len(failed_reads)} ({len(failed_reads)/len(common_reads)*100:.1f}%)")
    print()
    print(f"EXTEND could fix: {fixed_by_extend}/{len(failed_reads)} wrong reads")
    print()
    
    # Calculate expected NEW accuracy
    old_correct = len(common_reads) - len(failed_reads)
    new_correct = old_correct + fixed_by_extend
    old_accuracy = old_correct / len(common_reads) * 100
    new_accuracy = new_correct / len(common_reads) * 100
    improvement = new_accuracy - old_accuracy
    
    print(f"Expected improvement:")
    print(f"  OLD accuracy: {old_accuracy:.1f}%")
    print(f"  NEW accuracy (with EXTEND): {new_accuracy:.1f}%")
    print(f"  Improvement: +{improvement:.1f}%")
    print()
    
    if improvement > 20:
        print("🎉 EXTEND phase shows MAJOR improvement!")
    elif improvement > 10:
        print("✅ EXTEND phase shows good improvement")
    elif improvement > 0:
        print("✓ EXTEND phase helps somewhat")
    else:
        print("⚠️  EXTEND phase doesn't help")

if __name__ == '__main__':
    test_extend_on_failed_reads()
