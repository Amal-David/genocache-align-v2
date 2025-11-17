#!/usr/bin/env python3
"""
Robust GRCh38 encoding with checkpointing and error handling
"""
import os, sys, random, time, gc
import numpy as np
from Bio import SeqIO
import torch
import torch.nn as nn
import torch.nn.functional as F

# Config - UPDATED FOR FULL GENOME MODEL
REF_FA = "/home/nebius/genocache/GRCh38.fa"
ENCODER_PATH = "/home/nebius/genocache/nal_encoder_full_genome.pt"  # NEW TRAINED MODEL
OUT_VECTORS = "/home/nebius/genocache/grch38_full_vectors.npy"
OUT_POSITIONS = "/home/nebius/genocache/grch38_full_positions.npy"
OUT_CHROMOSOMES = "/home/nebius/genocache/grch38_full_chromosomes.npy"
CHECKPOINT_DIR = "/home/nebius/genocache/encoding_checkpoints_full"
SEED_LEN = 256
EMB_DIM = 256
STRIDE = 32
BATCH = 2048
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Create checkpoint directory
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

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

def save_checkpoint(chr_name, chr_idx, vectors, positions, chromosomes):
    """Save chromosome checkpoint"""
    checkpoint_file = f"{CHECKPOINT_DIR}/chr_{chr_idx}_{chr_name}.npz"
    try:
        np.savez_compressed(
            checkpoint_file,
            vectors=vectors,
            positions=np.array(positions, dtype=np.int64),
            chromosomes=np.array(chromosomes, dtype=np.int32)
        )
        print(f"    ✓ Checkpoint saved: {checkpoint_file}")
        return True
    except Exception as e:
        print(f"    ✗ Checkpoint save failed: {e}")
        return False

def load_checkpoints():
    """Load all existing checkpoints"""
    import glob
    checkpoint_files = sorted(glob.glob(f"{CHECKPOINT_DIR}/chr_*.npz"))
    
    if not checkpoint_files:
        return None, None, None, 0
    
    print(f"Found {len(checkpoint_files)} checkpoint files")
    
    all_vecs, all_pos, all_chr = [], [], []
    max_idx = 0
    
    for ckpt in checkpoint_files:
        data = np.load(ckpt)
        all_vecs.append(data['vectors'])
        all_pos.extend(data['positions'].tolist())
        all_chr.extend(data['chromosomes'].tolist())
        
        # Extract chr_idx from filename
        import re
        match = re.search(r'chr_(\d+)_', os.path.basename(ckpt))
        if match:
            chr_idx = int(match.group(1))
            max_idx = max(max_idx, chr_idx)
    
    print(f"Loaded {len(all_vecs)} chromosomes from checkpoints")
    return all_vecs, all_pos, all_chr, max_idx + 1

def main():
    print("=" * 70)
    print("Robust GRCh38 Encoding with Checkpointing")
    print("=" * 70)
    
    if not os.path.exists(REF_FA):
        print(f"ERROR: {REF_FA} not found")
        sys.exit(1)
    
    # Load model
    print(f"\n[1/3] Loading encoder: {ENCODER_PATH}")
    try:
        model = ImprovedNAL().to(DEVICE)
        model.load_state_dict(torch.load(ENCODER_PATH, map_location=DEVICE))
        model.eval()
        print(f"  ✓ Model loaded on {DEVICE}")
    except Exception as e:
        print(f"  ✗ Failed to load model: {e}")
        sys.exit(1)
    
    # Check for existing checkpoints
    print(f"\n[2/3] Checking for existing checkpoints...")
    all_vectors, all_positions, all_chromosomes, start_chr_idx = load_checkpoints()
    
    if all_vectors is None:
        all_vectors, all_positions, all_chromosomes = [], [], []
        start_chr_idx = 0
        print(f"  Starting fresh (no checkpoints found)")
    else:
        print(f"  ✓ Resuming from chromosome {start_chr_idx}")
    
    # Process genome
    print(f"\n[3/3] Processing chromosomes...")
    print(f"  Window size: {SEED_LEN}bp, Stride: {STRIDE}bp")
    
    chr_idx = 0
    total_bp = 0
    total_windows = 0
    start_time = time.time()
    
    # Track skipped/failed chromosomes
    skipped_chromosomes = []
    failed_chromosomes = []
    
    try:
        for record in SeqIO.parse(REF_FA, "fasta"):
            chr_name = record.id
            
            # Skip if already processed
            if chr_idx < start_chr_idx:
                chr_idx += 1
                continue
            
            seq = str(record.seq).upper()
            chr_len = len(seq)
            
            # Skip if too short or not main chromosome
            if chr_len < SEED_LEN:
                skipped_chromosomes.append(f"{chr_name} (too short: {chr_len} bp)")
                continue
            if not (chr_name.startswith('NC_') or chr_name.startswith('chr')):
                skipped_chromosomes.append(f"{chr_name} (not main chromosome)")
                continue
            
            print(f"\n  [{chr_idx}] Processing {chr_name}: {chr_len:,} bp")
            
            # Extract windows
            windows = []
            positions = []
            
            for i in range(0, chr_len - SEED_LEN + 1, STRIDE):
                win = seq[i:i+SEED_LEN]
                if win.count('N') < SEED_LEN * 0.2:
                    windows.append(win)
                    positions.append(i)
            
            if not windows:
                print(f"    Skipped (too many Ns)")
                skipped_chromosomes.append(f"{chr_name} (too many Ns)")
                chr_idx += 1
                continue
            
            print(f"    {len(windows):,} windows extracted")
            print(f"    Encoding...", end='', flush=True)
            
            # Encode in batches
            chr_vectors = []
            try:
                with torch.no_grad():
                    for i in range(0, len(windows), BATCH):
                        batch = windows[i:i+BATCH]
                        arr = np.stack([seq_to_onehot(s) for s in batch])
                        tensor = torch.from_numpy(arr).to(DEVICE)
                        vecs = model(tensor).cpu().numpy()
                        chr_vectors.append(vecs)
                        
                        # Clear GPU cache periodically
                        if i % (BATCH * 10) == 0:
                            torch.cuda.empty_cache() if torch.cuda.is_available() else None
                
                chr_vectors = np.vstack(chr_vectors)
                print(f" done ({chr_vectors.shape[0]:,} vectors)")
                
                # Save checkpoint immediately
                chr_positions = np.array(positions, dtype=np.int64)
                chr_chromosomes = np.array([chr_idx] * len(positions), dtype=np.int32)
                
                if save_checkpoint(chr_name, chr_idx, chr_vectors, chr_positions, chr_chromosomes):
                    # Add to global lists AFTER successful checkpoint
                    all_vectors.append(chr_vectors)
                    all_positions.extend(positions)
                    all_chromosomes.extend([chr_idx] * len(positions))
                    
                    chr_idx += 1
                    total_bp += chr_len
                    total_windows += len(windows)
                    
                    elapsed = time.time() - start_time
                    rate = total_bp / elapsed / 1e6
                    print(f"    Progress: {total_bp/1e9:.2f} Gbp @ {rate:.1f} Mbp/s")
                    
                    # Force garbage collection
                    del chr_vectors, chr_positions, chr_chromosomes, windows, positions
                    gc.collect()
                else:
                    print(f"    ✗ Failed to save checkpoint, skipping chromosome")
                    failed_chromosomes.append(f"{chr_name} (checkpoint save failed)")
                    chr_idx += 1
                    
            except Exception as e:
                print(f"\n    ✗ Error encoding {chr_name}: {e}")
                print(f"    Skipping to next chromosome...")
                failed_chromosomes.append(f"{chr_name} (encoding error: {str(e)[:50]})")
                chr_idx += 1
                continue
        
        # Final assembly from checkpoints
        print(f"\n[Final] Assembling final files from checkpoints...")
        
        import glob
        checkpoint_files = sorted(glob.glob(f"{CHECKPOINT_DIR}/chr_*.npz"))
        print(f"  Found {len(checkpoint_files)} checkpoint files")
        
        final_vectors = []
        final_positions = []
        final_chromosomes = []
        
        for i, ckpt in enumerate(checkpoint_files):
            print(f"  Loading {i+1}/{len(checkpoint_files)}: {os.path.basename(ckpt)}")
            data = np.load(ckpt)
            final_vectors.append(data['vectors'])
            final_positions.extend(data['positions'].tolist())
            final_chromosomes.extend(data['chromosomes'].tolist())
            
            # Free memory
            del data
            gc.collect()
        
        print(f"\n  Concatenating {len(final_vectors)} chromosome arrays...")
        vectors = np.vstack(final_vectors).astype('float32')
        positions = np.array(final_positions, dtype=np.int64)
        chromosomes = np.array(final_chromosomes, dtype=np.int32)
        
        print(f"  ✓ Total: {len(vectors):,} vectors")
        print(f"  ✓ Shape: {vectors.shape}")
        
        # Save final files in chunks to avoid memory issues
        print(f"\n  Saving final outputs...")
        
        print(f"    Saving {OUT_VECTORS}...")
        np.save(OUT_VECTORS, vectors)
        print(f"    ✓ Saved ({vectors.nbytes / 1024**3:.2f} GB)")
        
        print(f"    Saving {OUT_POSITIONS}...")
        np.save(OUT_POSITIONS, positions)
        print(f"    ✓ Saved ({positions.nbytes / 1024**2:.1f} MB)")
        
        print(f"    Saving {OUT_CHROMOSOMES}...")
        np.save(OUT_CHROMOSOMES, chromosomes)
        print(f"    ✓ Saved ({chromosomes.nbytes / 1024**2:.1f} MB)")
        
        elapsed = time.time() - start_time
        print(f"\n{'='*70}")
        print(f"Encoding complete in {elapsed/60:.1f} minutes")
        print(f"Average: {total_bp / elapsed / 1e6:.1f} Mbp/s")
        print(f"Total vectors: {len(vectors):,}")
        print(f"Total chromosomes: {chr_idx}")
        
        # Report skipped/failed chromosomes
        if skipped_chromosomes:
            print(f"\nSkipped chromosomes ({len(skipped_chromosomes)}):")
            for skip in skipped_chromosomes:
                print(f"  - {skip}")
        
        if failed_chromosomes:
            print(f"\n⚠️  Failed chromosomes ({len(failed_chromosomes)}):")
            for fail in failed_chromosomes:
                print(f"  ✗ {fail}")
            
            # Save failure log
            with open("encoding_failures.log", "w") as f:
                f.write("Failed Chromosomes:\n")
                for fail in failed_chromosomes:
                    f.write(f"{fail}\n")
            print(f"\n  Failure details saved to: encoding_failures.log")
        
        print("=" * 70)
        
    except KeyboardInterrupt:
        print(f"\n\n✗ Interrupted by user")
        print(f"Progress saved in {CHECKPOINT_DIR}/")
        print(f"Re-run this script to resume from checkpoint")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n✗ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        print(f"\nProgress saved in {CHECKPOINT_DIR}/")
        print(f"Re-run this script to resume from checkpoint")
        sys.exit(1)

if __name__ == "__main__":
    main()
