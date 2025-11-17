#!/usr/bin/env python3
"""
Autonomous GenoCache V4 Pipeline

This script runs the entire validation and scaling pipeline automatically:
1. Monitor chr22 training completion
2. Generate test reads
3. Run validation
4. If accuracy >80%, scale to full genome
5. If accuracy <80%, save debug info and wait

No human intervention needed if things go well!
"""

import subprocess
import time
import sys
from pathlib import Path
from datetime import datetime
import json


def log(msg):
    """Print with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")
    sys.stdout.flush()


def wait_for_training_completion(log_file: Path, max_wait_hours: float = 2.0):
    """Wait for chr22 training to complete"""
    log("=" * 80)
    log("STEP 1: Monitoring chr22 training")
    log("=" * 80)
    
    start_time = time.time()
    max_wait_seconds = max_wait_hours * 3600
    
    last_epoch = 0
    
    while True:
        elapsed = time.time() - start_time
        
        if elapsed > max_wait_seconds:
            log(f"❌ Training exceeded {max_wait_hours}h timeout")
            return False
        
        # Check if training is complete
        if log_file.exists():
            with open(log_file, 'r') as f:
                content = f.read()
                
                # Check for completion
                if "Training Complete!" in content:
                    log("✅ chr22 training completed!")
                    return True
                
                # Extract current epoch
                lines = content.split('\n')
                for line in reversed(lines[-100:]):
                    if line.startswith("Epoch "):
                        try:
                            epoch = int(line.split('/')[0].split()[-1])
                            if epoch > last_epoch:
                                last_epoch = epoch
                                log(f"  Training progress: Epoch {epoch}/50 ({epoch*2}%)")
                        except:
                            pass
                        break
        
        # Wait before checking again
        time.sleep(60)  # Check every minute


def generate_test_reads():
    """Generate chr22 test reads"""
    log("\n" + "=" * 80)
    log("STEP 2: Generating chr22 test reads")
    log("=" * 80)
    
    script = Path("/home/nebius/genocache/genocache-v4/validation/generate_chr22_test_reads.py")
    
    result = subprocess.run(
        f"cd /home/nebius/genocache/genocache-v4/validation && "
        f"source /home/nebius/genocache/.venv/bin/activate && "
        f"python3 {script}",
        shell=True,
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        log("✅ Test reads generated successfully")
        return True
    else:
        log(f"❌ Test read generation failed:")
        log(result.stderr)
        return False


def run_validation():
    """Run validation on chr22"""
    log("\n" + "=" * 80)
    log("STEP 3: Running chr22 validation")
    log("=" * 80)
    
    script = Path("/home/nebius/genocache/genocache-v4/validation/validate_chr22.py")
    
    result = subprocess.run(
        f"cd /home/nebius/genocache/genocache-v4/validation && "
        f"source /home/nebius/genocache/.venv/bin/activate && "
        f"python3 {script}",
        shell=True,
        capture_output=True,
        text=True
    )
    
    log(result.stdout)
    
    if result.returncode != 0:
        log(f"❌ Validation failed:")
        log(result.stderr)
        return None
    
    # Parse accuracy from output
    for line in result.stdout.split('\n'):
        if 'Accuracy @ ±1kb:' in line:
            try:
                accuracy = float(line.split(':')[1].strip().replace('%', ''))
                log(f"✅ Validation complete: {accuracy:.2f}% accuracy")
                return accuracy
            except:
                pass
    
    log("⚠️  Could not parse accuracy from output")
    return None


def create_full_genome_training_script():
    """Create full genome training script"""
    log("\n" + "=" * 80)
    log("STEP 4: Creating full genome training script")
    log("=" * 80)
    
    script_content = '''#!/usr/bin/env python3
"""
GenoCache V4 - Full Genome Training (Scale-up after chr22 success)

Training on all 25 chromosomes with proven configuration
"""

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from tqdm import tqdm
import time
from datetime import datetime

sys.path.append(str(Path(__file__).parent.parent))
from models.encoder import GenoCacheEncoder
sys.path.append(str(Path(__file__).parent))
from augmentation import DataAugmentation


class InfoNCELoss(nn.Module):
    def __init__(self, temperature=0.07):
        super().__init__()
        self.temperature = temperature
    
    def forward(self, anchor_emb, positive_emb):
        batch_size = anchor_emb.shape[0]
        anchor_emb = F.normalize(anchor_emb, p=2, dim=1)
        positive_emb = F.normalize(positive_emb, p=2, dim=1)
        similarity = torch.mm(anchor_emb, positive_emb.T) / self.temperature
        pos_sim = torch.diagonal(similarity)
        log_sum_exp = torch.logsumexp(similarity, dim=1)
        loss = -torch.mean(pos_sim - log_sum_exp)
        pos_sim_avg = torch.mean(pos_sim).item()
        mask = torch.eye(batch_size, device=similarity.device).bool()
        neg_similarities = similarity.masked_fill(mask, float('-inf'))
        neg_sim_avg = torch.mean(neg_similarities[~mask]).item()
        return loss, pos_sim_avg, neg_sim_avg


def load_genome(fasta_path: Path) -> dict:
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


def generate_training_batch(genome_dict, augmenter, batch_size=1024, device='cuda'):
    # Sample from random chromosome
    chr_names = list(genome_dict.keys())
    chr_weights = [len(seq) for seq in genome_dict.values()]
    chr_name = np.random.choice(chr_names, p=np.array(chr_weights)/sum(chr_weights))
    genome_seq = genome_dict[chr_name]
    
    anchors, positives = augmenter.create_training_pairs(genome_seq, batch_size=batch_size)
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def seq_to_tensor(seq):
        indices = [base_to_idx.get(b, 4) for b in seq]
        return torch.tensor(indices, dtype=torch.long)
    
    anchor_tensors = torch.stack([seq_to_tensor(seq) for seq in anchors])
    positive_tensors = torch.stack([seq_to_tensor(seq) for seq in positives])
    
    return anchor_tensors.to(device), positive_tensors.to(device)


def train_epoch(model, genome_dict, augmenter, criterion, optimizer, device, 
                batches_per_epoch=500, batch_size=1024):
    model.train()
    epoch_loss = 0
    epoch_pos_sim = 0
    epoch_neg_sim = 0
    
    pbar = tqdm(range(batches_per_epoch), desc="Training")
    for batch_idx in pbar:
        anchor_batch, positive_batch = generate_training_batch(
            genome_dict, augmenter, batch_size, device
        )
        
        anchor_emb = model(anchor_batch)
        positive_emb = model(positive_batch)
        loss, pos_sim, neg_sim = criterion(anchor_emb, positive_emb)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        epoch_loss += loss.item()
        epoch_pos_sim += pos_sim
        epoch_neg_sim += neg_sim
        
        separation = pos_sim - neg_sim
        pbar.set_postfix({
            'loss': f'{loss.item():.4f}',
            'sep': f'{separation:.4f}',
            'pos': f'{pos_sim:.4f}',
            'neg': f'{neg_sim:.4f}'
        })
    
    n_batches = batches_per_epoch
    return {
        'loss': epoch_loss / n_batches,
        'pos_sim': epoch_pos_sim / n_batches,
        'neg_sim': epoch_neg_sim / n_batches,
        'separation': (epoch_pos_sim - epoch_neg_sim) / n_batches
    }


def main():
    print("=" * 80)
    print("GenoCache V4 - Full Genome Training")
    print("=" * 80)
    print()
    
    GENOME_PATH = Path("/home/nebius/genocache/GRCh38.fa")
    CHECKPOINT_DIR = Path("/home/nebius/genocache/genocache-v4/models/checkpoints")
    LOG_DIR = Path("/home/nebius/genocache/genocache-v4/logs")
    
    # Increased for full genome
    BATCH_SIZE = 1024
    BATCHES_PER_EPOCH = 500  # 500 batches × 1024 = 512K examples per epoch
    NUM_EPOCHS = 30  # Less epochs needed with more data
    LEARNING_RATE = 1e-4
    TEMPERATURE = 0.07
    
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    LOG_DIR.mkdir(exist_ok=True)
    
    print("Configuration:")
    print(f"  Chromosomes: ALL (25)")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Batches per epoch: {BATCHES_PER_EPOCH}")
    print(f"  Examples per epoch: {BATCH_SIZE * BATCHES_PER_EPOCH:,}")
    print(f"  Epochs: {NUM_EPOCHS}")
    print(f"  Total examples: {BATCH_SIZE * BATCHES_PER_EPOCH * NUM_EPOCHS:,}")
    print()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
    print()
    
    # Load genome
    genome = load_genome(GENOME_PATH)
    print(f"✅ Loaded {len(genome)} chromosomes")
    total_bp = sum(len(seq) for seq in genome.values())
    print(f"  Total: {total_bp:,} bp")
    print()
    
    # Create augmenter
    augmenter = DataAugmentation(seed_len=512)
    print("✅ Augmenter ready\\n")
    
    # Load chr22 model as warm start
    print("Loading chr22 model as warm start...")
    chr22_checkpoints = list(CHECKPOINT_DIR.glob("chr22_best_*.pt"))
    if not chr22_checkpoints:
        chr22_checkpoints = list(CHECKPOINT_DIR.glob("chr22_epoch*.pt"))
    
    model = GenoCacheEncoder(
        emb_dim=128,
        seed_len=512,
        vocab_size=5,
        hidden_dims=[64, 128, 256],
        num_attention_layers=2,
        dropout=0.1
    ).to(device)
    
    if chr22_checkpoints:
        checkpoint_path = sorted(chr22_checkpoints)[-1]
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        print(f"✅ Loaded chr22 checkpoint: {checkpoint_path.name}")
    else:
        print("⚠️  No chr22 checkpoint found, training from scratch")
    
    num_params = sum(p.numel() for p in model.parameters())
    print(f"✅ Model ready: {num_params:,} parameters\\n")
    
    criterion = InfoNCELoss(temperature=TEMPERATURE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    
    print("Starting training...\\n")
    
    best_separation = 0
    start_time = time.time()
    
    log_file = LOG_DIR / f"training_fullgenome_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    with open(log_file, 'w') as f:
        f.write("epoch,loss,pos_sim,neg_sim,separation\\n")
        
        for epoch in range(1, NUM_EPOCHS + 1):
            print(f"\\n{'='*80}")
            print(f"Epoch {epoch}/{NUM_EPOCHS}")
            print(f"{'='*80}")
            
            train_metrics = train_epoch(
                model, genome, augmenter, criterion, optimizer, device,
                batches_per_epoch=BATCHES_PER_EPOCH,
                batch_size=BATCH_SIZE
            )
            
            print(f"\\nEpoch {epoch} Summary:")
            print(f"  Loss: {train_metrics['loss']:.4f}")
            print(f"  Pos sim: {train_metrics['pos_sim']:.4f}")
            print(f"  Neg sim: {train_metrics['neg_sim']:.4f}")
            print(f"  Separation: {train_metrics['separation']:.4f}")
            
            f.write(f"{epoch},{train_metrics['loss']:.6f},{train_metrics['pos_sim']:.6f},"
                   f"{train_metrics['neg_sim']:.6f},{train_metrics['separation']:.6f}\\n")
            f.flush()
            
            # Save checkpoint
            if train_metrics['separation'] > best_separation:
                best_separation = train_metrics['separation']
                checkpoint_path = CHECKPOINT_DIR / f"fullgenome_best_sep{train_metrics['separation']:.4f}_epoch{epoch}.pt"
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'train_metrics': train_metrics,
                }, checkpoint_path)
                print(f"\\n✅ Saved best checkpoint: {checkpoint_path.name}")
            
            if epoch % 5 == 0:
                checkpoint_path = CHECKPOINT_DIR / f"fullgenome_epoch{epoch}.pt"
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'train_metrics': train_metrics,
                }, checkpoint_path)
                print(f"💾 Saved checkpoint: {checkpoint_path.name}")
            
            elapsed = time.time() - start_time
            avg_epoch_time = elapsed / epoch
            remaining = avg_epoch_time * (NUM_EPOCHS - epoch)
            print(f"\\nTime: {elapsed/3600:.1f}h elapsed, {remaining/3600:.1f}h remaining")
    
    total_time = time.time() - start_time
    print("\\n" + "="*80)
    print("Full Genome Training Complete!")
    print("="*80)
    print(f"Total time: {total_time/3600:.2f} hours")
    print(f"Best separation: {best_separation:.4f}")
    print(f"Log file: {log_file}")


if __name__ == "__main__":
    main()
'''
    
    script_path = Path("/home/nebius/genocache/genocache-v4/training/train_fullgenome.py")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    script_path.chmod(0o755)
    log(f"✅ Created: {script_path}")
    return script_path


def start_full_genome_training(script_path: Path):
    """Start full genome training"""
    log("\n" + "=" * 80)
    log("STEP 5: Starting full genome training")
    log("=" * 80)
    
    log_file = Path("/home/nebius/genocache/genocache-v4/logs/train_fullgenome.out")
    pid_file = Path("/home/nebius/genocache/genocache-v4/logs/train_fullgenome.pid")
    
    cmd = (
        f"cd /home/nebius/genocache/genocache-v4/training && "
        f"source /home/nebius/genocache/.venv/bin/activate && "
        f"nohup python3 -u {script_path} > {log_file} 2>&1 & "
        f"echo $! > {pid_file} && echo $!"
    )
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
    
    if result.returncode == 0:
        pid = result.stdout.strip()
        log(f"✅ Full genome training started! PID: {pid}")
        log(f"  Log: {log_file}")
        log(f"  Estimated time: 10-15 hours")
        return True
    else:
        log(f"❌ Failed to start training:")
        log(result.stderr)
        return False


def save_summary(accuracy: float, decision: str):
    """Save summary for user"""
    summary = {
        'timestamp': datetime.now().isoformat(),
        'chr22_accuracy': accuracy,
        'decision': decision,
        'target': 80.0,
        'status': 'success' if accuracy >= 80 else 'needs_attention'
    }
    
    summary_file = Path("/home/nebius/genocache/genocache-v4/PIPELINE_SUMMARY.json")
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    log(f"\n✅ Summary saved: {summary_file}")


def main():
    log("=" * 80)
    log("GenoCache V4 - Autonomous Pipeline")
    log("=" * 80)
    log("Starting autonomous pipeline...")
    log("This will run steps 1-5 automatically based on chr22 results")
    log("=" * 80)
    
    # Step 1: Wait for chr22 training
    chr22_log = Path("/home/nebius/genocache/genocache-v4/logs/train_chr22.out")
    if not wait_for_training_completion(chr22_log):
        log("\n❌ PIPELINE FAILED: Training timeout or error")
        return 1
    
    # Step 2: Generate test reads
    if not generate_test_reads():
        log("\n❌ PIPELINE FAILED: Test read generation")
        return 1
    
    # Step 3: Run validation
    accuracy = run_validation()
    if accuracy is None:
        log("\n❌ PIPELINE FAILED: Validation error")
        return 1
    
    # Step 4: Decision point
    if accuracy >= 80:
        log(f"\n🎉 SUCCESS! Accuracy {accuracy:.2f}% >= 80%")
        log("Proceeding to full genome training automatically...")
        
        script_path = create_full_genome_training_script()
        if start_full_genome_training(script_path):
            save_summary(accuracy, "scaled_to_full_genome")
            log("\n" + "=" * 80)
            log("PIPELINE COMPLETE - Full genome training in progress")
            log("=" * 80)
            log(f"chr22 accuracy: {accuracy:.2f}%")
            log("Full genome training: RUNNING")
            log("Estimated completion: 10-15 hours")
            log("\nYou can sleep! Everything is running autonomously.")
            return 0
        else:
            log("\n❌ Failed to start full genome training")
            return 1
    
    elif accuracy >= 60:
        log(f"\n⚠️  MARGINAL: Accuracy {accuracy:.2f}% (60-80%)")
        log("Needs analysis - waiting for human review")
        save_summary(accuracy, "marginal_needs_analysis")
        return 2
    
    else:
        log(f"\n❌ FAILED: Accuracy {accuracy:.2f}% < 60%")
        log("Needs debugging - waiting for human review")
        save_summary(accuracy, "failed_needs_debug")
        return 3


if __name__ == "__main__":
    sys.exit(main())
