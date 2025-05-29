#!/bin/sh

# Default settings
BUILD_TYPE="Release"
CXX_FLAGS=""
C_FLAGS=""
USE_OMP="ON"

# Parse first argument
if [ "$1" = "debug" ]; then
    BUILD_TYPE="Debug"
    CXX_FLAGS="-ggdb3 -fsanitize=undefined -fsanitize=address"
    C_FLAGS="-ggdb3 -fsanitize=undefined -fsanitize=address"
elif [ "$1" = "valgrind" ]; then
    BUILD_TYPE="Debug"
    CXX_FLAGS="-ggdb3"
    C_FLAGS="-ggdb3"
elif [ "$1" = "release" ] || [ -z "$1" ]; then
    BUILD_TYPE="Release"
    CXX_FLAGS=""
    C_FLAGS=""
fi

echo "Building executable with configuration: ${BUILD_TYPE}"
echo "OpenMP enabled: ${USE_OMP}"


# Cannot be done here because this will only change the environment variable locally ...
#if [ "${USE_OMP}" = "ON" ]; then
#	export OMP_NUM_THREADS="$(nproc --all)"
#fi

mkdir -p build
cd build 

cmake \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" \
    -DCMAKE_C_FLAGS="${C_FLAGS}" \
    -DUSE_OMP="${USE_OMP}" \
    ../src_tlb

make -j"$(nproc)"
