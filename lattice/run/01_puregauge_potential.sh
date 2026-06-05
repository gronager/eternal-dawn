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

# --- 1. generate the pure-gauge ensemble --------------------------------------------------
# NB: Grid's Test_hmc_WilsonGauge fixes beta and trajectory count in its source and writes
# NERSC checkpoints named ckpoint_lat.<n> in the CURRENT directory. Run it from $OUT, and set
# beta/trajectories via its real flags -- check '$HMC --help'. The bare run already thermalises.
( cd "$OUT" && "$HMC" --grid "$GRIDSPEC" --accelerator-threads 8 2>&1 | tee hmc.log )

# --- 2. measure W(R,T) on every saved config ----------------------------------------------
: > "$OUT/wloops_raw.dat"        # columns: R  T  W   (one block per config)
shopt -s nullglob
cfgs=("$OUT"/ckpoint_lat.*)
[ "${#cfgs[@]}" -gt 0 ] || { echo "no ckpoint_lat.* in $OUT -- did the HMC write configs?"; exit 1; }
for cfg in "${cfgs[@]}"; do
  [[ "$cfg" == *.gz ]] && continue
  # keep only the R T W data lines; Grid prints init chatter to stdout too (it starts with
  # letters: "MPI...", "Grid...", "AcceleratorCudaInit...", "local...", "SharedMemory...").
  "$MEAS" --grid "$GRIDSPEC" --config "$cfg" --rmax $((L/2)) --tmax $((T/2)) \
      --accelerator-threads 8 | grep -E '^[0-9]' >> "$OUT/wloops_raw.dat"
done

# --- 3. average over configs, extract V(R), fit the string tension -------------------------
echo "== analyze: string tension =="
python3 - "$OUT/wloops_raw.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
# robust: keep only clean R T W triples (in case any Grid stdout chatter slipped through)
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
# average W over configs for each (R,T)
keys = sorted(set(map(tuple, raw[:,:2].astype(int))))
Ru, Tu, Wu = [], [], []
for (r,t) in keys:
    m = (R==r)&(T==t)
    Ru.append(r); Tu.append(t); Wu.append(W[m].mean())
Rr, Vr = lat.effective_potential(np.array(Ru), np.array(Tu), np.array(Wu))
good = np.isfinite(Vr)
fit = lat.static_potential_cornell(Rr[good], Vr[good])
print(f"  sigma  = {fit['sigma']:.5f} a^-2   (>0 => confinement / area law; ~0 => screening)")
print(f"  alpha  = {fit['alpha']:.4f}   r0(Sommer) = {fit['r0_sommer']:.3f} a")
PY
