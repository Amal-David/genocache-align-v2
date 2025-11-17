"""
Build FAISS indexes from reference embeddings
Supports: Flat (validation), IVFFlat, IVFPQ, OPQ+IVFPQ (production)
"""

import faiss
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import argparse


def build_flat_index(vectors, output_path):
    """Build uncompressed flat index for validation"""
    print("\nBuilding Flat index (exact search)...")
    
    d = vectors.shape[1]
    index = faiss.IndexFlatIP(d)  # Inner product (for normalized vectors = cosine similarity)
    index.add(vectors)
    
    faiss.write_index(index, str(output_path))
    
    print(f"✅ Built Flat index")
    print(f"   Vectors: {index.ntotal:,}")
    print(f"   Size: {output_path.stat().st_size / 1e6:.1f} MB")
    
    return index


def build_ivf_index(vectors, output_path, nlist=4096, nprobe=32):
    """Build IVF partitioned index"""
    print(f"\nBuilding IVF index (nlist={nlist})...")
    
    d = vectors.shape[1]
    quantizer = faiss.IndexFlatIP(d)
    index = faiss.IndexIVFFlat(quantizer, d, nlist)
    
    # Train the index
    print("  Training quantizer...")
    index.train(vectors)
    
    # Add vectors
    print("  Adding vectors...")
    index.add(vectors)
    
    # Set search parameters
    index.nprobe = nprobe
    
    faiss.write_index(index, str(output_path))
    
    print(f"✅ Built IVF index")
    print(f"   Vectors: {index.ntotal:,}")
    print(f"   Clusters: {nlist}")
    print(f"   nprobe: {nprobe}")
    print(f"   Size: {output_path.stat().st_size / 1e6:.1f} MB")
    
    return index


def build_ivfpq_index(vectors, output_path, nlist=4096, m=32, nbits=8, nprobe=32):
    """Build IVF + Product Quantization index (compressed)"""
    print(f"\nBuilding IVFPQ index (nlist={nlist}, m={m}, nbits={nbits})...")
    
    d = vectors.shape[1]
    quantizer = faiss.IndexFlatIP(d)
    index = faiss.IndexIVFPQ(quantizer, d, nlist, m, nbits)
    
    # Train the index
    print("  Training quantizer and PQ...")
    index.train(vectors)
    
    # Add vectors
    print("  Adding vectors...")
    index.add(vectors)
    
    # Set search parameters
    index.nprobe = nprobe
    
    faiss.write_index(index, str(output_path))
    
    compression_ratio = (d * 4) / ((m * nbits) / 8)
    
    print(f"✅ Built IVFPQ index")
    print(f"   Vectors: {index.ntotal:,}")
    print(f"   Clusters: {nlist}")
    print(f"   PQ subvectors: {m}")
    print(f"   Bits per subvector: {nbits}")
    print(f"   Compression: {compression_ratio:.1f}x")
    print(f"   nprobe: {nprobe}")
    print(f"   Size: {output_path.stat().st_size / 1e6:.1f} MB")
    
    return index


def build_opq_ivfpq_index(vectors, output_path, nlist=4096, m=32, nbits=8, nprobe=32):
    """Build OPQ + IVF + PQ index (best compression)"""
    print(f"\nBuilding OPQ+IVFPQ index (nlist={nlist}, m={m}, nbits={nbits})...")
    
    d = vectors.shape[1]
    
    # OPQ requires d to be divisible by m
    if d % m != 0:
        print(f"  Warning: d={d} not divisible by m={m}, adjusting...")
        m_adj = d // (d // m)
        print(f"  Using m={m_adj} instead")
        m = m_adj
    
    # Build OPQ + IVFPQ index
    quantizer = faiss.IndexFlatIP(d)
    index_base = faiss.IndexIVFPQ(quantizer, d, nlist, m, nbits)
    
    # Add OPQ preprocessing
    index = faiss.IndexPreTransform(
        faiss.OPQMatrix(d, m),
        index_base
    )
    
    # Train the index
    print("  Training OPQ, quantizer, and PQ...")
    index.train(vectors)
    
    # Add vectors
    print("  Adding vectors...")
    index.add(vectors)
    
    # Set search parameters
    index_base.nprobe = nprobe
    
    faiss.write_index(index, str(output_path))
    
    compression_ratio = (d * 4) / ((m * nbits) / 8)
    
    print(f"✅ Built OPQ+IVFPQ index")
    print(f"   Vectors: {index.ntotal:,}")
    print(f"   Clusters: {nlist}")
    print(f"   PQ subvectors: {m}")
    print(f"   Bits per subvector: {nbits}")
    print(f"   Compression: {compression_ratio:.1f}x")
    print(f"   nprobe: {nprobe}")
    print(f"   Size: {output_path.stat().st_size / 1e6:.1f} MB")
    
    return index


def build_indexes(vectors_path, manifest_path, output_dir, index_types=['flat', 'ivfpq']):
    """
    Build multiple FAISS indexes from reference embeddings
    
    Args:
        vectors_path: Path to ref_vectors.npy
        manifest_path: Path to ref_manifest.json
        output_dir: Where to save indexes
        index_types: List of index types to build
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*80}")
    print("BUILDING FAISS INDEXES")
    print(f"{'='*80}")
    print(f"Vectors: {vectors_path}")
    print(f"Output: {output_dir}")
    print(f"Index types: {index_types}")
    print(f"{'='*80}\n")
    
    # Load vectors
    print("Loading reference vectors...")
    vectors = np.load(vectors_path).astype('float32')
    print(f"✅ Loaded {len(vectors):,} vectors ({vectors.shape[1]}-D)")
    
    # Ensure vectors are normalized for cosine similarity
    print("\nNormalizing vectors...")
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    vectors = vectors / (norms + 1e-8)
    print(f"✅ Normalized (mean norm: {np.linalg.norm(vectors, axis=1).mean():.4f})")
    
    # Load manifest
    with open(manifest_path) as f:
        manifest = json.load(f)
    
    # Build indexes
    indexes_built = {}
    
    if 'flat' in index_types:
        flat_path = output_dir / 'index_flat.faiss'
        build_flat_index(vectors, flat_path)
        indexes_built['flat'] = str(flat_path)
    
    if 'ivf' in index_types:
        ivf_path = output_dir / 'index_ivf.faiss'
        build_ivf_index(vectors, ivf_path, nlist=4096, nprobe=32)
        indexes_built['ivf'] = str(ivf_path)
    
    if 'ivfpq' in index_types:
        ivfpq_path = output_dir / 'index_ivfpq.faiss'
        build_ivfpq_index(vectors, ivfpq_path, nlist=4096, m=32, nbits=8, nprobe=32)
        indexes_built['ivfpq'] = str(ivfpq_path)
    
    if 'opq' in index_types:
        opq_path = output_dir / 'index_opq_ivfpq.faiss'
        build_opq_ivfpq_index(vectors, opq_path, nlist=4096, m=32, nbits=8, nprobe=32)
        indexes_built['opq'] = str(opq_path)
    
    # Save index manifest
    index_manifest = {
        'created_at': datetime.now().isoformat(),
        'source': {
            'vectors_path': str(vectors_path),
            'manifest_path': str(manifest_path),
            'checkpoint': manifest['checkpoint'],
            'reference': manifest['reference']
        },
        'indexes': indexes_built,
        'num_vectors': len(vectors),
        'embedding_dim': vectors.shape[1]
    }
    
    index_manifest_path = output_dir / 'indexes_manifest.json'
    with open(index_manifest_path, 'w') as f:
        json.dump(index_manifest, f, indent=2)
    
    print(f"\n✅ Saved index manifest: {index_manifest_path}")
    
    print(f"\n{'='*80}")
    print("INDEX BUILDING COMPLETE")
    print(f"{'='*80}")
    print(f"Indexes built: {len(indexes_built)}")
    for idx_type, path in indexes_built.items():
        size = Path(path).stat().st_size / 1e6
        print(f"  {idx_type}: {size:.1f} MB")
    print(f"{'='*80}\n")
    
    return indexes_built


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build FAISS indexes')
    
    parser.add_argument('--vectors', type=str, required=True,
                       help='Path to reference vectors (.npy)')
    parser.add_argument('--manifest', type=str, required=True,
                       help='Path to reference manifest (.json)')
    parser.add_argument('--output-dir', type=str, default='genocache_v2/indexes',
                       help='Output directory for indexes')
    parser.add_argument('--index-types', type=str, nargs='+',
                       default=['flat', 'ivfpq'],
                       choices=['flat', 'ivf', 'ivfpq', 'opq'],
                       help='Types of indexes to build')
    
    args = parser.parse_args()
    
    build_indexes(
        vectors_path=args.vectors,
        manifest_path=args.manifest,
        output_dir=args.output_dir,
        index_types=args.index_types
    )
