#!/usr/bin/env bash
# Build the Grid lattice library on a GH200 (Grace/aarch64 + Hopper sm_90), single-GPU.
# Grid (https://github.com/paboyle/Grid) ships its own GPU Dirac solvers; QUDA is optional
# for the high-performance multigrid solver (see the note at the bottom).
#
# This is a STARTING template -- adapt paths/module names to your box. It assumes a CUDA
# toolkit (>=12.x for Hopper) and a recent C++ host compiler are available.
set -euo pipefail

PREFIX="${PREFIX:-$HOME/ed-lattice}"          # install prefix for deps + Grid
NPROC="${NPROC:-$(nproc)}"
CUDA_ARCH="${CUDA_ARCH:-sm_90}"               # Hopper
mkdir -p "$PREFIX/src" && cd "$PREFIX/src"

echo "== 1. dependencies (GMP, MPFR, FFTW, c-lime, HDF5) =="
# Most are available via apt/spack/conda; install to $PREFIX. Example with the system pkgs:
#   sudo apt-get install -y libgmp-dev libmpfr-dev libfftw3-dev libhdf5-dev
# c-lime (for ILDG/SciDAC I/O):
if [ ! -d c-lime ]; then
  git clone https://github.com/usqcd-software/c-lime.git
  ( cd c-lime && ./autogen.sh && ./configure --prefix="$PREFIX" && make -j"$NPROC" && make install )
fi

echo "== 2. clone Grid =="
[ -d Grid ] || git clone https://github.com/paboyle/Grid.git
cd Grid && ./bootstrap.sh && mkdir -p build && cd build

echo "== 3. configure Grid for single-GPU Hopper (no MPI) =="
# Key flags: CUDA accelerator, GPU SIMD, no inter-node comms (single GH200),
# nvcc as the device compiler with the Hopper gencode.
../configure \
  --prefix="$PREFIX" \
  --enable-comms=none \
  --enable-simd=GPU \
  --enable-accelerator=cuda \
  --enable-gen-simd-width=64 \
  --enable-unified=yes \
  CXX=nvcc \
  CXXFLAGS="-ccbin g++ -gencode arch=compute_${CUDA_ARCH#sm_},code=${CUDA_ARCH} -std=c++17 -O3 --expt-relaxed-constexpr --expt-extended-lambda -I$PREFIX/include" \
  LDFLAGS="-L$PREFIX/lib" \
  --with-lime="$PREFIX" --with-fftw=/usr --with-gmp=/usr --with-mpfr=/usr --with-hdf5=/usr

echo "== 4. build =="
make -j"$NPROC"
make install

echo "== done. Grid installed under $PREFIX =="
echo "Smoke test:  $PREFIX/src/Grid/build/tests/Test_dwf_hdcr --grid 8.8.8.8 --accelerator-threads 8"
echo
echo "OPTIONAL high-performance solver: build QUDA with"
echo "  cmake -DQUDA_GPU_ARCH=${CUDA_ARCH} -DQUDA_DIRAC_WILSON=ON -DQUDA_DIRAC_CLOVER=ON \\"
echo "        -DQUDA_MULTIGRID=ON -DQUDA_DOWNLOAD_USQCD=ON .. && make -j"
echo "and reconfigure Grid with --enable-quda=<quda-prefix> to offload the solves."
