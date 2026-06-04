#!/usr/bin/env bash
# Target 1: string tension sigma from the static potential (pure SU(3) gauge, single GH200).
# Generate a pure-gauge ensemble (Wilson action), measure Wilson loops -> V(r), fit Cornell.
# Cost: hours. Validates the stack and confirms the area law (confinement) before any fermions.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # sets GRID, require_exe, note_params
OUT="${OUT:-out/puregauge}"; mkdir -p "$OUT"
require_exe "$GRID/tests/hmc/Test_hmc_WilsonGauge"
require_exe "$GRID/tests/core/Test_WilsonLoops"
note_params "$GRID/tests/hmc/Test_hmc_WilsonGauge"

L=24 ; T=48            # 24^3 x 48 -- comfortably single-GPU
BETA=6.0               # SU(3) Wilson beta (a ~ 0.09 fm here; adjust for your scale)
NCFG=200 ; NSEP=10     # measure every 10 trajectories

# --- generate the pure-gauge ensemble (heatbath/HMC) ---
"$GRID/tests/hmc/Test_hmc_WilsonGauge" \
    --grid $L.$L.$L.$T --beta $BETA \
    --Trajectories $((NCFG*NSEP)) --Thermalizations 100 \
    --Checkpoint "$OUT/cfg" --accelerator-threads 8  2>&1 | tee "$OUT/hmc.log"

# --- measure the static potential from smeared Wilson loops on each saved config ---
# (Test_WilsonLoops / a smearing+loop measurement; writes r, V(r) pairs to potential.dat)
for cfg in "$OUT"/cfg.*; do
    "$GRID/tests/core/Test_WilsonLoops" --config "$cfg" --grid $L.$L.$L.$T \
        --rmax $((L/2)) >> "$OUT/potential_raw.dat"
done
# average over configs -> r, V, V_err  (one row per r)
python3 - "$OUT/potential_raw.dat" "$OUT/potential.dat" <<'PY'
import sys, numpy as np
raw = np.loadtxt(sys.argv[1])          # columns: r  V  (per config)
r = np.unique(raw[:,0])
rows = [[rr, raw[raw[:,0]==rr,1].mean(), raw[raw[:,0]==rr,1].std(ddof=1)/np.sqrt((raw[:,0]==rr).sum())] for rr in r]
np.savetxt(sys.argv[2], rows, header="r  V  V_err")
PY

echo "== analyze: string tension =="
python3 - "$OUT/potential.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
r, V, Verr = np.loadtxt(sys.argv[1], unpack=True)
fit = lat.static_potential_cornell(r, V, Verr)
print(f"  sigma  = {fit['sigma']:.5f} a^-2   (>0 => confinement / area law)")
print(f"  alpha  = {fit['alpha']:.4f}   r0(Sommer) = {fit['r0_sommer']:.3f} a")
PY
