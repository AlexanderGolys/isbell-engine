#!/bin/bash

# Build script for Linux (Arch)
# Usage: ./rebuild.sh [SCENE_NAME]

set -e

# Default scene
SCENE_NAME="${1:-impulsiveInteractions}"

# Check if premake5 is available
if ! command -v premake5 &> /dev/null; then
    echo "premake5 not found. Installing via pacman..."
    sudo pacman -S premake --noconfirm
fi

# Clean previous build if it exists
if [ -d "build" ]; then
    echo "Cleaning previous build..."
    rm -rf build
fi

# Generate CMake project with premake5
echo "Generating CMake project for scene: $SCENE_NAME"
premake5 --scene="$SCENE_NAME" cmake

# Build with CMake
echo "Building project..."
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel$(nproc)

echo "Build complete. Binary located at: build/$SCENE_NAME/bin/build-x64/"
