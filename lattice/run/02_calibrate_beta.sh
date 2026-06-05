#!/usr/bin/env bash
# Calibrate the bare coupling against the plaquette (Eternal Dawn lattice).
# The input "beta" fed to the gauge action need not equal the standard 6/g^2 -- the equilibrium
# is set by the action, and the only convention-free anchor is a MEASURED observable. The SU(3)
# Wilson reference is <P> = 0.5669 at beta = 5.6. We scan several input betas (one per stream, in
# parallel under MPS), let each equilibrate, measure the equilibrium plaquette, and interpolate
# the input beta whose plaquette hits the target. That calibrated beta is what 00/01 should use.
#
#   ./lattice/run/02_calibrate_beta.sh                 # default scan around 5.6..6.2
#   BETAS="5.8 5.9 6.0 6.1" TARGET=0.5669 ./lattice/run/02_calibrate_beta.sh
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, start_mps/stop_mps, run_pool
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/calib}"; mkdir -p "$OUT"
GEN="$HERE/../src/generate_gauge"
MEAS="$HERE/../src/measure_potential"
require_exe "$GEN"; require_exe "$MEAS"

L="${L:-16}" ; T="${T:-32}" ; GRIDSPEC="$L.$L.$L.$T"
BETAS="${BETAS:-5.7 5.9 6.0 6.2}"   # input betas to scan (plaquette is monotone in beta)
NTRAJ="${NTRAJ:-80}"                # plaquette is UV-fast: ~40 traj to equilibrate, measure the rest
THERM="${THERM:-40}"                # discard below this trajectory (plaquette equilibration)
TARGET="${TARGET:-0.5669}"          # SU(3) Wilson reference plaquette at standard beta=5.6

echo "calibrating beta against the plaquette (target <P>=$TARGET): scanning betas [$BETAS]"
# generate only the betas that don't already have enough configs (a re-run reuses prior streams)
gen=""; i=0
for b in $BETAS; do
  i=$((i+1)); d="$OUT/beta$b"
  have=$(ls -v "$d"/ckpoint_lat.* 2>/dev/null | grep -vcE '\.gz$|rng' || true)
  if [ "$have" -ge "$NTRAJ" ]; then
    echo "  beta=$b: $have configs already present -- skipping generation"
  else
    gen+=$(printf 'mkdir -p %q && cd %q && %q --grid %q --beta %q --seed %d --StartingType HotStart --Trajectories %d --accelerator-threads 8 > hmc.log 2>&1\n' \
      "$d" "$d" "$GEN" "$GRIDSPEC" "$b" "$((i*1000 + 1))" "$NTRAJ")$'\n'
  fi
done
if [ -n "$gen" ]; then
  start_mps
  printf '%s' "$gen" | run_pool "$(echo $BETAS | wc -w)"   # all missing betas side by side under MPS
  stop_mps
fi

# report what generation produced before measuring (so a failed/slow generation is obvious, not silent)
echo "== generated configs per beta =="
for b in $BETAS; do
  d="$OUT/beta$b"
  cnt=$(ls -v "$d"/ckpoint_lat.* 2>/dev/null | grep -vcE '\.gz$|rng' || true)
  echo "  beta=$b  $cnt configs in $d  (see $d/hmc.log if 0)"
done

# measure the equilibrium plaquette of each beta (average over thermalised configs).
# Defensive against set -e/pipefail: extraction can return empty without aborting the run.
plaq_of() {  # plaq_of <config>  -> our own avgPlaquette line, or empty
  "$MEAS" --grid "$GRIDSPEC" --config "$1" --rmax 1 --tmax 1 --accelerator-threads 8 2>/dev/null \
    | grep -m1 '^# config' | sed -n 's/.*plaquette //p' | awk '{print $1}' || true
}
echo "== equilibrium plaquette per beta =="
: > "$OUT/plaq_vs_beta.dat"
for b in $BETAS; do
  d="$OUT/beta$b"
  vals=()
  for c in $(ls -v "$d"/ckpoint_lat.* 2>/dev/null | grep -vE '\.gz$|rng' || true); do
    n="${c##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
    (( n >= THERM )) || continue
    p="$(plaq_of "$c")"
    [ -n "$p" ] && vals+=("$p")
  done
  if [ "${#vals[@]}" -gt 0 ]; then
    pbar=$(printf '%s\n' "${vals[@]}" | awk '{s+=$1; n++} END{printf "%.6f", s/n}')
    echo "  beta=$b  <P>=$pbar  (over ${#vals[@]} configs)"
    echo "$b $pbar" >> "$OUT/plaq_vs_beta.dat"
  else
    echo "  beta=$b  -- no thermalised configs >= traj $THERM (raise NTRAJ, or check $d/hmc.log)"
  fi
done

echo "== calibrated beta =="
TARGET="$TARGET" python3 - "$OUT/plaq_vs_beta.dat" <<'PY'
import os, sys, numpy as np
from cartasis_sims import lattice as lat
target = float(os.environ["TARGET"])
data = np.loadtxt(sys.argv[1], ndmin=2)
if data.size == 0 or data.shape[0] < 2:
    sys.exit("need >=2 betas with measured plaquettes -- raise NTRAJ or add BETAS")
betas, plaqs = data[:,0], data[:,1]
out = lat.beta_from_plaquette(betas, plaqs, target=target)
how = "interpolated" if out["interpolated"] else "EXTRAPOLATED (scan a wider BETAS range to confirm)"
print(f"  measured: " + "  ".join(f"b={b:.3f}:<P>={p:.4f}" for b, p in zip(betas, plaqs)))
print(f"  -> input beta giving <P>={target}:  beta* = {out['beta']:.3f}   [{how}]")
print(f"  use this beta* in 00_generate.sh / 01_puregauge_potential.sh for true beta=5.6 physics")
PY
