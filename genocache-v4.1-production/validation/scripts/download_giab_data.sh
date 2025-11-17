#!/bin/bash
#
# Download GIAB HG002 ONT reads for validation
#
# Options:
# 1. Direct download from ONT Open Data (AWS S3)
# 2. From NCBI SRA
# 3. Using existing local data

set -e

echo "============================================================"
echo "GIAB HG002 Data Download for Validation"
echo "============================================================"
echo ""

# Check for existing data
EXISTING_DATA="/home/nebius/genocache/giab_hg002_chr22_ont.fastq.gz"

if [ -f "$EXISTING_DATA" ]; then
    echo "✅ Found existing GIAB data: $EXISTING_DATA"
    
    # Check if it's empty
    if [ -s "$EXISTING_DATA" ]; then
        echo "   File size: $(du -h $EXISTING_DATA | cut -f1)"
        echo ""
        echo "Use this file? (y/n)"
        read -r response
        if [[ "$response" == "y" ]]; then
            # Extract subset for validation
            echo "Extracting 100 reads for validation..."
            zcat "$EXISTING_DATA" | head -400 > ../data/giab_hg002_chr22_100reads.fastq
            echo "✅ Created validation subset: ../data/giab_hg002_chr22_100reads.fastq"
            exit 0
        fi
    else
        echo "⚠️  File exists but is empty"
    fi
fi

echo "Need to download GIAB HG002 data"
echo ""
echo "Options:"
echo "  1. AWS S3 (ONT Open Data) - RECOMMENDED"
echo "  2. Download using wget"
echo "  3. Skip and use existing data"
echo ""
echo "Choice (1-3):"
read -r choice

case $choice in
    1)
        echo "Downloading from AWS S3..."
        echo ""
        
        # Check if AWS CLI is available
        if command -v aws &> /dev/null; then
            echo "Using AWS CLI..."
            aws s3 ls s3://ont-open-data/giab_2025.01/ --no-sign-request | head -20
            
            echo ""
            echo "Note: Full dataset is very large (>100GB)"
            echo "We'll download a smaller subset for validation"
            echo ""
            echo "Proceed? (y/n)"
            read -r proceed
            
            if [[ "$proceed" == "y" ]]; then
                # This would download the full data - adjust as needed
                # aws s3 cp s3://ont-open-data/giab_2025.01/HG002/ ../data/ --recursive --no-sign-request
                echo "⚠️  Please manually download a subset from AWS S3"
                echo "   or use option 2 for smaller downloads"
            fi
        else
            echo "❌ AWS CLI not installed"
            echo "   Install with: pip install awscli"
        fi
        ;;
    
    2)
        echo "Downloading using wget..."
        echo ""
        echo "Note: This is a placeholder - actual URL needs to be specified"
        echo "Real GIAB data URLs:"
        echo "  - https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/"
        echo "  - Check for HG002 ONT reads"
        ;;
    
    3)
        echo "Skipping download"
        echo ""
        echo "Using placeholder..."
        echo "⚠️  No real data available for validation"
        ;;
    
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "============================================================"
echo "Download Status"
echo "============================================================"
if [ -f "../data/giab_hg002_chr22_100reads.fastq" ]; then
    echo "✅ Validation data ready"
    wc -l ../data/giab_hg002_chr22_100reads.fastq
else
    echo "⚠️  No validation data available yet"
    echo "   Please download manually or use existing local data"
fi
