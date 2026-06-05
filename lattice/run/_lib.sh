# Shared helpers for the lattice run scripts. Source this at the top of each.
# Sets GRID (the Grid build dir) and provides require_exe to fail loudly + helpfully.
GRID="${GRID:-$HOME/ed-lattice/src/Grid/build}"

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
