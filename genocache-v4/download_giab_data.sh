#!/bin/bash
#
# Download GIAB HG002 ONT reads for real data testing
# Source: Genome in a Bottle Consortium
#

set -e

DATA_DIR="/home/nebius/genocache/genocache_data/giab"
mkdir -p "$DATA_DIR"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║        Downloading GIAB HG002 ONT Reads for Testing                ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

# HG002 ONT ultra-long reads from GIAB
# Using chr22 subset for quick testing
GIAB_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/"

# Download small subset for testing
echo "Downloading HG002 ONT reads (subset)..."
echo "Note: Full dataset is ~100GB, downloading small test set"
echo ""

# For now, let's use a smaller ONT dataset from AWS
# AWS Public Datasets - ONT Human
AWS_ONT_URL="s3://ont-open-data/gm24385_2020.09/"

echo "Alternative: Using AWS ONT open data..."
echo "Checking if AWS CLI is available..."

if command -v aws &> /dev/null; then
    echo "✅ AWS CLI found"
    echo ""
    echo "Downloading chr22 reads subset (this may take 10-30 min)..."
    
    # Download just chr22 reads if available
    # For testing, we'll download a small FASTQ subset
    aws s3 cp --no-sign-request \
        s3://ont-open-data/gm24385_2020.09/ \
        "$DATA_DIR/" \
        --recursive \
        --exclude "*" \
        --include "*chr22*" \
        --include "*fastq.gz" \
        2>/dev/null || echo "Direct chr22 download not available"
    
else
    echo "⚠️  AWS CLI not found"
    echo ""
    echo "Alternative: Using wget to download from NCBI..."
    
    # Download a small subset using wget
    # Using Rel15 data which has chr-specific BAMs
    NCBI_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/"
    
    echo "Downloading small test dataset from NCBI..."
    cd "$DATA_DIR"
    
    # Download a small FASTQ file for testing (~1-2 GB)
    wget -c "${NCBI_URL}140818_D00360_0047_AHA66FADXX/Project_RM8398_L8537_CAGATC/Sample_2A/2A_CAGATC_L004_R1_001.fastq.gz" \
        -O hg002_test_r1.fastq.gz 2>&1 | tail -20
    
    if [ $? -eq 0 ]; then
        echo "✅ Download complete!"
        echo ""
        ls -lh "$DATA_DIR"/*.fastq.gz
    else
        echo "⚠️  Download failed, trying alternative..."
    fi
fi

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    ALTERNATIVE: Generate Test Reads                 ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "If downloads fail, we can use our existing synthetic reads or"
echo "generate more realistic test reads with error profiles."
echo ""
echo "Data directory: $DATA_DIR"
echo ""
