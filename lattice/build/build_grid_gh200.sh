#!/usr/bin/env bash
# Build Grid + its test executables on a GH200 (Grace/aarch64 + Hopper sm_90), single-GPU,
# AND install every external dependency it needs. Idempotent: re-running skips finished steps.
#
# Fixes the "Test_hmc_WilsonGauge: No such file or directory" trap: Grid's `make` builds only
# the LIBRARY and the top-level tests; the run scripts use tests in sub-directories
# (tests/hmc, tests/core, ...), which are built only by `make tests`. We run it.
set -euo pipefail

# Absolute path to the repo's lattice/ dir -- captured NOW, before any `cd`, because later
# steps cd into the Grid build tree and $0 may be relative (that broke step 6 before).
LATTICE_DIR="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"

PREFIX="${PREFIX:-$HOME/ed-lattice}"
NPROC="${NPROC:-$(nproc)}"
CUDA_ARCH="${CUDA_ARCH:-sm_90}"                 # Hopper (GH200)
GRID_BRANCH="${GRID_BRANCH:-develop}"
mkdir -p "$PREFIX/src"

# --------------------------------------------------------------------------
echo "== 1. external dependencies =="
# Grid needs: autotools, GMP, MPFR, FFTW3, HDF5, OpenSSL, zlib, and c-lime (SciDAC/ILDG I/O).
APT_PKGS="build-essential automake autoconf libtool pkg-config git \
          libgmp-dev libmpfr-dev libfftw3-dev libhdf5-dev libssl-dev zlib1g-dev"
if command -v apt-get >/dev/null 2>&1; then
  SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo"
  $SUDO apt-get update -y
  $SUDO apt-get install -y $APT_PKGS
else
  echo "!! apt-get not found. Install these yourself before continuing:"
  echo "   $APT_PKGS"
  echo "   (on conda:  conda install -c conda-forge gmp mpfr fftw hdf5 openssl autoconf automake libtool)"
fi

# c-lime is not packaged; build from source into $PREFIX.
if [ ! -f "$PREFIX/lib/liblime.a" ]; then
  echo "-- building c-lime --"
  cd "$PREFIX/src"
  [ -d c-lime ] || git clone https://github.com/usqcd-software/c-lime.git
  cd c-lime && ./autogen.sh && ./configure --prefix="$PREFIX" \
    && make -j"$NPROC" && make install
fi

# --------------------------------------------------------------------------
echo "== 2. clone Grid ($GRID_BRANCH) =="
cd "$PREFIX/src"
[ -d Grid ] || git clone https://github.com/paboyle/Grid.git
cd Grid && git checkout "$GRID_BRANCH" && git pull --ff-only || true
[ -f configure ] || ./bootstrap.sh
mkdir -p build && cd build

# --------------------------------------------------------------------------
echo "== 3. configure for single-GPU Hopper (no MPI) =="
if [ ! -f Makefile ]; then
  ../configure \
    --prefix="$PREFIX" \
    --enable-comms=none \
    --enable-simd=GPU \
    --enable-accelerator=cuda \
    --enable-gen-simd-width=64 \
    --enable-unified=yes \
    --with-lime="$PREFIX" \
    --with-gmp=/usr --with-mpfr=/usr --with-fftw=/usr --with-hdf5=/usr \
    CXX=nvcc \
    CXXFLAGS="-ccbin g++ -gencode arch=compute_${CUDA_ARCH#sm_},code=${CUDA_ARCH} -std=c++17 -O3 --expt-relaxed-constexpr --expt-extended-lambda -I$PREFIX/include" \
    LDFLAGS="-L$PREFIX/lib"
fi

# --------------------------------------------------------------------------
echo "== 4. build the library, install it, then build ONLY the executables we need =="
make -j"$NPROC"          # the Grid library (+ top-level tests)
make install             # headers + libGrid under $PREFIX

# Grid's full `make tests` is all-or-nothing: ONE sibling test that does not compile on a given
# (GPU) config aborts the whole build -- which is exactly the tests/core error you saw. We don't
# need the suite, only five executables, so build each one independently; a broken sibling can no
# longer block us, and a wrong name (Grid versions differ) just prints a warning, not a halt.
build_one() {   # build_one <subdir> <exe>
  echo "-- building $1/$2 --"
  if make -C "$1" -j"$NPROC" "$2"; then
    echo "   ok: $1/$2"
  else
    echo "!! $1/$2 did NOT build -- its name may differ in this Grid version (see step 5),"
    echo "   or it has a genuine compile error; the OTHER targets below are unaffected."
  fi
}
build_one tests/hmc      Test_hmc_WilsonGauge
build_one tests/hmc      Test_hmc_WilsonFermionGauge
build_one tests/core     Test_WilsonLoops
build_one tests/smearing Test_WilsonFlow
build_one tests/solver   Test_ModeNumber

# --------------------------------------------------------------------------
echo "== 5. verify the executables the run scripts need actually exist =="
GRID_BUILD="$PWD"
need=(
  "tests/hmc/Test_hmc_WilsonGauge"
  "tests/hmc/Test_hmc_WilsonFermionGauge"
  "tests/core/Test_WilsonLoops"
  "tests/smearing/Test_WilsonFlow"
)
missing=0
for t in "${need[@]}"; do
  if [ -x "$GRID_BUILD/$t" ]; then echo "  ok   $t"; else echo "  MISS $t"; missing=1; fi
done
echo
echo "All built HMC/measurement tests (use these exact names in the run scripts):"
find "$GRID_BUILD/tests" -maxdepth 2 -type f -perm -u+x -name 'Test_*' 2>/dev/null \
  | sed "s|$GRID_BUILD/||" | sort
echo
echo "GRID build dir: $GRID_BUILD"
echo "Point the run scripts at it with:  export GRID=$GRID_BUILD"
if [ "$missing" -ne 0 ]; then
  echo "!! Some expected executables are missing -- Grid's test names vary by version."
  echo "   Pick the matching ones from the list above and edit lattice/run/*.sh accordingly."
fi

# --------------------------------------------------------------------------
echo "== 6. build our custom measurement programs (static potential, ...) =="
if [ -d "$LATTICE_DIR/src" ]; then
  make -C "$LATTICE_DIR/src" PREFIX="$PREFIX" GRID_CONFIG="$PREFIX/bin/grid-config" \
    && echo "   ok: measurement programs in $LATTICE_DIR/src" \
    || echo "!! measurement programs did not build -- run 'make -C lattice/src' and paste the error"
else
  echo "(lattice/src not found beside the build script; skipping custom programs)"
fi
echo "== done =="
