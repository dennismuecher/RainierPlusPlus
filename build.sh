#!/bin/bash
set -e

echo "Building RAINIER 2.0..."
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc 2>/dev/null || echo 4)
echo "Build complete! Executable: $(pwd)/rainier"
