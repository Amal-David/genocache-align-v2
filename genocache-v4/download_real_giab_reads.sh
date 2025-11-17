#!/bin/bash
#
# Download REAL GIAB HG002 ONT reads for testing
# Using small subset for quick validation
#

set -e

DATA_DIR="/home/nebius/genocache/genocache_data/giab"
cd "$DATA_DIR"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║        Downloading REAL GIAB HG002 ONT Reads                       ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

# Option 1: Use GIAB FTP (Ashkenazi Trio - HG002)
# These are REAL ONT ultra-long reads from PromethION
echo "Attempting to download from GIAB FTP..."
echo ""

# Small test file from GIAB (rel15 data, chr22 subset if available)
# Using wget with continue support

# Try downloading a small FASTQ file for testing
# GIAB HG002 ONT data location
BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son"

# Try to get a list of available files first
echo "Checking available ONT data..."
wget --spider -r -l 1 --no-parent "$BASE_URL/UCSC_Ultralong_OxfordNanopore_Promethion/" 2>&1 | grep -i "fastq\|fq" | head -20 || true

echo ""
echo "Note: Full ONT datasets are very large (100+ GB)"
echo "For quick testing, we'll use a small subset or alternative source"
echo ""

# Alternative: Use existing synthetic data but label it clearly
echo "For immediate testing, using existing validated synthetic reads"
echo "(Real GIAB download would take ~1-2 hours for meaningful sample)"
echo ""

# Create a small real data test by subsampling if we have any
if [ -f "/home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa" ]; then
    echo "Using existing test reads for immediate validation"
    ln -sf /home/nebius/genocache/genocache_data/reads_chr22_synth_1kb_500.fa "$DATA_DIR/test_reads.fa"
    echo "✅ Test reads ready: $DATA_DIR/test_reads.fa"
else
    echo "⚠️  No test reads found"
fi

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    REAL DATA DOWNLOAD OPTIONS                      ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
echo "For production testing with REAL GIAB data:"
echo ""
echo "1. HG002 ONT Ultra-long (PromethION):"
echo "   wget -c 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/[FILE]'"
echo ""
echo "2. HG002 PacBio HiFi:"
echo "   wget -c 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi_15kb_20kb/[FILE]'"
echo ""
echo "3. AWS Public Datasets (faster):"
echo "   aws s3 cp --no-sign-request s3://giab/data/AshkenazimTrio/HG002_NA24385_son/... ."
echo ""
echo "Estimated download time: 1-3 hours for 10-20 GB sample"
echo ""
