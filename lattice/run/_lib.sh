# Shared helpers for the lattice run scripts. Source this at the top of each.
# Sets GRID (the Grid build dir) and provides require_exe to fail loudly + helpfully.
GRID="${GRID:-$HOME/ed-lattice/src/Grid/build}"

# Make the repo's Python package (cartasis_sims) importable in the analysis steps without
# requiring a pip install on the box -- the run scripts call `from cartasis_sims import lattice`.
_LATTICE_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PYTHONPATH="$_LATTICE_REPO/sims/src${PYTHONPATH:+:$PYTHONPATH}"

require_exe() {
  # require_exe <path-to-executable>
  local exe="$1"
  if [ ! -x "$exe" ]; then
    echo "ERROR: required Grid executable not found:" >&2
    echo "   $exe" >&2
    echo >&2
    echo "Build Grid AND its tests first (the build script now runs 'make tests'):" >&2
    echo "   lattice/build/build_grid_gh200.sh" >&2
    echo >&2
    local dir; dir="$(dirname "$exe")"
    if [ -d "$dir" ]; then
      echo "Executables that DO exist in $dir :" >&2
      find "$dir" -maxdepth 1 -type f -perm -u+x -printf '   %f\n' 2>/dev/null >&2 \
        || ls "$dir" >&2
      echo "(if this is a Grid test dir, names vary by version; if it is lattice/src," >&2
      echo " run the build -- step 6 compiles our measurement programs.)" >&2
    else
      echo "(directory $dir does not exist -- the build did not complete.)" >&2
    fi
    exit 1
  fi
}

note_params() {
  # Grid's stock HMC tests sometimes fix the action parameters (e.g. beta) in the source
  # rather than on the command line. If a run ignores a flag below, check the binary's help
  # and, if needed, edit the .cc and rebuild, or supply Grid's HMC parameter/XML input.
  echo "[note] verify flags with '$1 --help'; some Grid HMC tests read parameters from source"
  echo "       or an XML input rather than CLI -- adjust if a parameter below is ignored."
}

# ---- packing many small jobs onto one big GPU (96 GB GH200) -------------------------------
# The lattices here use ~1-2 GB each, so the GPU is memory-rich but compute-underutilised by a
# single small job. MPS lets independent CUDA processes truly share the SMs (not just
# time-slice), so N parallel measurements run ~N x faster until the GPU saturates.
start_mps() {
  command -v nvidia-cuda-mps-control >/dev/null 2>&1 || { echo "[mps] not available; skipping"; return 0; }
  export CUDA_MPS_PIPE_DIRECTORY="${CUDA_MPS_PIPE_DIRECTORY:-/tmp/nvidia-mps}"
  export CUDA_MPS_LOG_DIRECTORY="${CUDA_MPS_LOG_DIRECTORY:-/tmp/nvidia-mps-log}"
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  if nvidia-cuda-mps-control -d 2>/dev/null; then
    echo "[mps] daemon started ($CUDA_MPS_PIPE_DIRECTORY)"
  else
    echo "[mps] daemon already running -- reusing it"   # already up: fine, not an error
  fi
  return 0                                               # never fail the caller (set -e safe)
}
stop_mps() {
  command -v nvidia-cuda-mps-control >/dev/null 2>&1 || return 0
  echo quit | nvidia-cuda-mps-control 2>/dev/null || true
  echo "[mps] daemon stopped"
  return 0
}

run_pool() {
  # run_pool <NPAR> -- read commands (one per line) from stdin, run NPAR at a time.
  local npar="$1"
  while IFS= read -r cmd; do
    bash -c "$cmd" &
    while [ "$(jobs -rp | wc -l)" -ge "$npar" ]; do wait -n; done
  done
  wait
}
