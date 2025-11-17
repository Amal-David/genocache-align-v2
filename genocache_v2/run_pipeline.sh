#!/bin/bash
# GenoCache-Align V2: Complete Pipeline
# Runs: Training → Encoding → Index Building

set -e
export PYTHONPATH="$(pwd)/scripts:${PYTHONPATH:-}"

# Configuration
REFERENCE_FASTA="references/GCF_000001405.40_GRCh38.p14_chr22.fna"
CHROM="NC_000022.11"
OUTPUT_BASE="genocache_v2"
NUM_GPUS=1

# Training parameters
EPOCHS=50
BATCH_SIZE=256
LEARNING_RATE=1e-3

# Encoding parameters
SEED_LEN=512
STRIDE=32
ENCODE_BATCH_SIZE=512

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║          GenoCache-Align V2: Complete Pipeline                ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"

# Activate environment
echo -e "\n${YELLOW}[1/5] Activating environment...${NC}"
source ~/miniconda/etc/profile.d/conda.sh
conda activate neuralign
echo -e "${GREEN}✓ Environment activated${NC}"

# Step 1: Prepare training data
echo -e "\n${YELLOW}[2/5] Preparing training data...${NC}"
python3 << EOF
from genomic_utils import create_training_data

create_training_data(
    reference_path='${REFERENCE_FASTA}',
    chrom='${CHROM}',
    seed_len=${SEED_LEN},
    stride=${STRIDE},
    output_dir='${OUTPUT_BASE}/data/training'
)
EOF
echo -e "${GREEN}✓ Training data ready${NC}"

# Step 2: Train model
echo -e "\n${YELLOW}[3/5] Training ImprovedCNN (${NUM_GPUS} GPUs, ${EPOCHS} epochs)...${NC}"
echo -e "${BLUE}This will take 1-2 hours...${NC}"

torchrun --nproc_per_node=${NUM_GPUS} train_fast_cnn.py \
    --data-path ${OUTPUT_BASE}/data/training/seeds_${CHROM}.npz \
    --output-dir ${OUTPUT_BASE}/models/cnn_fast \
    --epochs ${EPOCHS} \
    --batch-size ${BATCH_SIZE} \
    --learning-rate ${LEARNING_RATE}

echo -e "${GREEN}✓ Training complete${NC}"

# Step 3: Encode reference
echo -e "\n${YELLOW}[4/5] Encoding reference genome...${NC}"
python encode_reference.py \
    --checkpoint ${OUTPUT_BASE}/models/cnn_fast/checkpoints/improved_cnn_best.pt \
    --fasta ${REFERENCE_FASTA} \
    --chrom ${CHROM} \
    --output-dir ${OUTPUT_BASE}/data/reference \
    --seed-len ${SEED_LEN} \
    --stride ${STRIDE} \
    --batch-size ${ENCODE_BATCH_SIZE}

echo -e "${GREEN}✓ Reference encoded${NC}"

# Step 4: Build indexes
echo -e "\n${YELLOW}[5/5] Building FAISS indexes...${NC}"
python build_indexes.py \
    --vectors ${OUTPUT_BASE}/data/reference/ref_vectors_latest.npy \
    --manifest ${OUTPUT_BASE}/data/reference/ref_manifest_latest.json \
    --output-dir ${OUTPUT_BASE}/indexes/compressed \
    --index-types flat ivfpq

echo -e "${GREEN}✓ Indexes built${NC}"

# Summary
echo -e "\n${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                    PIPELINE COMPLETE!                         ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"

echo -e "\n${GREEN}Output files:${NC}"
echo "  Checkpoint:  ${OUTPUT_BASE}/models/cnn_fast/checkpoints/improved_cnn_best.pt"
echo "  Embeddings:  ${OUTPUT_BASE}/data/reference/ref_vectors_latest.npy"
echo "  Flat Index:  ${OUTPUT_BASE}/indexes/compressed/index_flat.faiss"
echo "  IVFPQ Index: ${OUTPUT_BASE}/indexes/compressed/index_ivfpq.faiss"

echo -e "\n${YELLOW}Next steps:${NC}"
echo "  1. Validate recall with test queries"
echo "  2. Tune FAISS nprobe parameter"
echo "  3. Scale to whole genome"
echo "  4. Add Hyena rescue tier"

echo -e "\n${BLUE}See GUIDE.md for validation examples!${NC}"
