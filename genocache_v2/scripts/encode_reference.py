"""
Encode reference genome into embeddings
Generates: ref_vectors.npy, ref_positions.npy, ref_manifest.json
"""

import torch
import numpy as np
from pathlib import Path
import json
from datetime import datetime
from tqdm import tqdm
import argparse

from improved_cnn import ImprovedCNN
from genomic_utils import ReferenceGenome, one_hot_encode


def encode_reference(checkpoint_path, fasta_path, chrom, output_dir, 
                     seed_len=512, stride=32, batch_size=512):
    """
    Encode reference genome into embeddings
    
    Args:
        checkpoint_path: Path to trained model checkpoint
        fasta_path: Path to reference FASTA
        chrom: Chromosome to encode
        output_dir: Where to save outputs
        seed_len: Seed length
        stride: Distance between seeds
        batch_size: Encoding batch size
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print(f"\n{'='*80}")
    print("ENCODING REFERENCE GENOME")
    print(f"{'='*80}")
    print(f"Checkpoint: {checkpoint_path}")
    print(f"Reference: {fasta_path}")
    print(f"Chromosome: {chrom}")
    print(f"Device: {device}")
    print(f"{'='*80}\n")
    
    # Load model
    print("Loading model...")
    model, metadata = ImprovedCNN.load_checkpoint(checkpoint_path, device=device)
    model.eval()
    
    emb_dim = metadata.get('emb_dim', model.out_dim)
    
    # Load reference
    ref = ReferenceGenome(fasta_path)
    
    # Extract seeds
    print(f"\nExtracting seeds from {chrom}...")
    seeds_data = ref.extract_seeds(chrom, seed_len, stride, both_strands=False)
    
    seeds = [s[0] for s in seeds_data]
    positions = [s[1] for s in seeds_data]
    
    print(f"✅ Extracted {len(seeds):,} seeds")
    print(f"   Genomic span: {min(positions):,} - {max(positions) + seed_len:,} bp")
    
    # Encode in batches
    print("\nEncoding seeds...")
    all_embeddings = []
    
    with torch.no_grad():
        for i in tqdm(range(0, len(seeds), batch_size)):
            batch_seeds = seeds[i:i + batch_size]
            
            # One-hot encode
            batch_tensors = torch.stack([one_hot_encode(s) for s in batch_seeds])
            batch_tensors = batch_tensors.to(device)
            
            # Encode
            embeddings = model(batch_tensors)
            all_embeddings.append(embeddings.cpu().numpy())
    
    # Concatenate all embeddings
    all_embeddings = np.vstack(all_embeddings)
    positions_array = np.array(positions, dtype=np.int64)
    
    print(f"\n✅ Encoded {len(all_embeddings):,} seeds")
    print(f"   Shape: {all_embeddings.shape}")
    print(f"   Dtype: {all_embeddings.dtype}")
    
    # Generate timestamp for versioning
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save embeddings
    vectors_path = output_dir / f'ref_vectors_{timestamp}.npy'
    np.save(vectors_path, all_embeddings)
    print(f"\n✅ Saved embeddings: {vectors_path}")
    
    # Save positions
    positions_path = output_dir / f'ref_positions_{timestamp}.npy'
    np.save(positions_path, positions_array)
    print(f"✅ Saved positions: {positions_path}")
    
    # Create manifest for provenance tracking
    manifest = {
        'created_at': datetime.now().isoformat(),
        'checkpoint': {
            'path': str(checkpoint_path),
            'timestamp': metadata.get('timestamp', 'unknown'),
            'architecture': metadata.get('architecture', 'ImprovedCNN')
        },
        'reference': {
            'fasta_path': str(fasta_path),
            'chromosome': chrom,
            'seed_length': seed_len,
            'stride': stride
        },
        'embeddings': {
            'path': str(vectors_path),
            'shape': list(all_embeddings.shape),
            'dtype': str(all_embeddings.dtype),
            'num_vectors': len(all_embeddings)
        },
        'positions': {
            'path': str(positions_path),
            'min': int(min(positions)),
            'max': int(max(positions)),
            'genomic_span': f"{min(positions)}-{max(positions) + seed_len}"
        }
    }
    
    manifest_path = output_dir / f'ref_manifest_{timestamp}.json'
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"✅ Saved manifest: {manifest_path}")
    
    # Also save as "latest" symlinks (for convenience)
    latest_vectors = output_dir / 'ref_vectors_latest.npy'
    latest_positions = output_dir / 'ref_positions_latest.npy'
    latest_manifest = output_dir / 'ref_manifest_latest.json'
    
    for src, dst in [(vectors_path, latest_vectors),
                     (positions_path, latest_positions),
                     (manifest_path, latest_manifest)]:
        if dst.exists():
            dst.unlink()
        dst.symlink_to(src.name)
    
    print(f"\n✅ Created 'latest' symlinks for easy access")
    
    print(f"\n{'='*80}")
    print("ENCODING COMPLETE")
    print(f"{'='*80}")
    print(f"Total vectors: {len(all_embeddings):,}")
    print(f"Embedding dim: {emb_dim}")
    print(f"File size: {vectors_path.stat().st_size / 1e6:.1f} MB")
    print(f"{'='*80}\n")
    
    return vectors_path, positions_path, manifest_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Encode reference genome')
    
    parser.add_argument('--checkpoint', type=str, required=True,
                       help='Path to trained model checkpoint')
    parser.add_argument('--fasta', type=str, required=True,
                       help='Path to reference FASTA file')
    parser.add_argument('--chrom', type=str, required=True,
                       help='Chromosome to encode (e.g., NC_000022.11)')
    parser.add_argument('--output-dir', type=str, default='genocache_v2/data/reference',
                       help='Output directory')
    parser.add_argument('--seed-len', type=int, default=512,
                       help='Seed length')
    parser.add_argument('--stride', type=int, default=32,
                       help='Distance between seeds')
    parser.add_argument('--batch-size', type=int, default=512,
                       help='Encoding batch size')
    
    args = parser.parse_args()
    
    encode_reference(
        checkpoint_path=args.checkpoint,
        fasta_path=args.fasta,
        chrom=args.chrom,
        output_dir=args.output_dir,
        seed_len=args.seed_len,
        stride=args.stride,
        batch_size=args.batch_size
    )
