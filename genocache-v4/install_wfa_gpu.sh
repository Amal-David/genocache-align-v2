#!/bin/bash
# Install WFA-GPU for production alignment

set -e

echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                    Installing WFA-GPU Library                              ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check CUDA availability
if ! command -v nvcc &> /dev/null; then
    echo "⚠️  CUDA not found. Checking nvidia-smi..."
    nvidia-smi || echo "No GPU detected"
fi

# Clone WFA-GPU repository
echo "Cloning WFA-GPU..."
if [ -d "WFA-GPU" ]; then
    echo "WFA-GPU directory already exists, pulling latest..."
    cd WFA-GPU && git pull && cd ..
else
    git clone https://github.com/quim0/WFA-GPU.git
fi

cd WFA-GPU

# Build library
echo ""
echo "Building WFA-GPU..."
make clean || true
make -j$(nproc)

echo ""
echo "✅ WFA-GPU build complete!"
echo ""
echo "Library location: $(pwd)/lib"
echo "Include location: $(pwd)/include"

cd ..

echo ""
echo "Next: Create Python bindings..."

