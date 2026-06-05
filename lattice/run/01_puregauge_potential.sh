#!/usr/bin/env bash
# Target 1: string tension sigma from the static potential (pure SU(3) gauge, single GH200).
# (1) generate a pure-gauge ensemble with Grid's HMC, (2) measure timelike Wilson loops W(R,T)
# with our measure_potential program, (3) extract V(R) and fit the Cornell form for sigma.
# Cost: hours. Validates the stack and confirms the area law (confinement) before any fermions.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # sets GRID, require_exe, note_params
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/puregauge}"; mkdir -p "$OUT"
HMC="$GRID/tests/hmc/Test_hmc_WilsonGauge"
MEAS="$HERE/../src/measure_potential"               # our custom program (built by the build script)
require_exe "$HMC"
require_exe "$MEAS"
note_params "$HMC"

L=16 ; T=32            # 16^3 x 32 -- a quick first ensemble (raise once it runs)
GRIDSPEC="$L.$L.$L.$T"
NTRAJ="${NTRAJ:-400}"  # total trajectories (Test_hmc_WilsonGauge: beta=5.6 hardcoded, saveInterval=1)
THERM="${THERM:-150}"  # discard configs below this trajectory (thermalisation)
STRIDE="${STRIDE:-20}" # then measure every STRIDE-th config (decorrelation)
NPAR="${NPAR:-8}"      # parallel measurement jobs on the GPU (tiny lattices -> pack the 96 GB)

# --- 1. generate the pure-gauge ensemble (skip if configs already exist) -------------------
# If you ran ./lattice/run/00_generate.sh first (recommended -- K independent streams in
# parallel), the configs are already in $OUT/stream*/ and this step is skipped. Otherwise it
# falls back to a single serial chain with the stock Test_hmc_WilsonGauge (beta=5.6 hardcoded).
shopt -s nullglob
if ls "$OUT"/ckpoint_lat.* "$OUT"/stream*/ckpoint_lat.* >/dev/null 2>&1; then
  echo "using existing configs in $OUT (run/00 streams or a previous run)"
else
  echo "no configs found -- generating a single serial chain (or run run/00 first for streams)"
  ( cd "$OUT" && "$HMC" --grid "$GRIDSPEC" --Trajectories "$NTRAJ" --accelerator-threads 8 \
      2>&1 | tee hmc.log )
fi

# --- 2. measure W(R,T) on THERMALISED, DECORRELATED configs (parallel on the GPU) ----------
sel=()
for cfg in "$OUT"/ckpoint_lat.* "$OUT"/stream*/ckpoint_lat.*; do
  [[ "$cfg" == *.gz ]] && continue
  n="${cfg##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
  (( n >= THERM )) || continue                 # thermalisation cut
  (( n % STRIDE == 0 )) || continue            # decorrelation stride
  sel+=("$cfg")
done
[ "${#sel[@]}" -gt 0 ] || { echo "no thermalised configs (n>=$THERM, stride $STRIDE) -- raise NTRAJ"; exit 1; }
echo "measuring ${#sel[@]} configs, $NPAR in parallel on the GPU"

start_mps                                      # true SM-sharing for the small jobs (no-op if absent)
# one self-contained measurement command per config (data lines only), run NPAR at a time
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --rmax %d --tmax %d --accelerator-threads 8 | grep -E "^[0-9]" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" $((L/2)) $((T/2)) "$cfg.wl"
done | run_pool "$NPAR"
stop_mps

cat "${sel[@]/%/.wl}" > "$OUT/wloops_raw.dat"  # combine the per-config results
echo "measured ${#sel[@]} configs"

# --- 3. average over configs, extract V(R), fit the string tension -------------------------
echo "== analyze: string tension =="
python3 - "$OUT/wloops_raw.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = []
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) == 3:
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            pass
raw = np.array(rows)
if raw.size == 0:
    sys.exit("no numeric R T W rows found in wloops_raw.dat")
R, T, W = raw[:,0], raw[:,1], raw[:,2]
Ru, Tu, Wu = [], [], []
for (r,t) in sorted(set(map(tuple, raw[:,:2].astype(int)))):
    m = (R==r)&(T==t)
    Ru.append(r); Tu.append(t); Wu.append(W[m].mean())
# robust V(R): slope of -ln W vs T over the small-T signal window (no smearing yet)
Rr, Vr = lat.potential_from_loops_fit(np.array(Ru), np.array(Tu), np.array(Wu), tmin=1, tmax=4)
good = np.isfinite(Vr) & (Vr > 0)
print(f"  V(R): R={Rr[good]}  V={np.round(Vr[good],4)}   (should rise with R)")
if good.sum() >= 3:
    fit = lat.static_potential_cornell(Rr[good], Vr[good])
    print(f"  sigma  = {fit['sigma']:.5f} a^-2   (>0 => confinement / area law; ~0 => screening)")
    print(f"  alpha  = {fit['alpha']:.4f}   r0(Sommer) = {fit['r0_sommer']:.3f} a")
else:
    print("  too few clean V(R) points to fit -- needs smearing and/or more configs")
PY
