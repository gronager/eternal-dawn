#!/usr/bin/env bash
# Target 2: deconfinement T_c from the Polyakov-loop susceptibility (pure gauge, single GH200).
# This is the lattice version of "the condensate melts" (genesis Stage 3): scan beta on a
# fixed-Nt lattice, find the susceptibility peak. T_c = 1/(Nt a(beta_c)); set a via target 3.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # sets GRID, require_exe, note_params
OUT="${OUT:-out/finiteT}"; mkdir -p "$OUT"
require_exe "$GRID/tests/hmc/Test_hmc_WilsonGauge"
note_params "$GRID/tests/hmc/Test_hmc_WilsonGauge"

NS=16 ; NT=6           # 16^3 x 6 : finite-T, Nt=6 sets the temperature axis
BETAS=$(seq 5.65 0.02 5.95)   # scan across the transition
NCFG=400

: > "$OUT/polyakov.dat"   # columns: beta  |L|  chi_L
for BETA in $BETAS; do
  "$GRID/tests/hmc/Test_hmc_WilsonGauge" \
      --grid $NS.$NS.$NS.$NT --beta $BETA \
      --Trajectories $NCFG --Thermalizations 100 \
      --PolyakovMeasure "$OUT/poly_$BETA.dat" --accelerator-threads 8 \
      2>&1 | tee -a "$OUT/hmc.log"
  python3 - "$BETA" "$OUT/poly_$BETA.dat" "$OUT/polyakov.dat" <<'PY'
import sys, numpy as np
beta=float(sys.argv[1]); L=np.loadtxt(sys.argv[2])         # per-config |L|
absL=np.abs(L); chi=absL.size*(np.mean(absL**2)-np.mean(absL)**2)
open(sys.argv[3],'a').write(f"{beta} {absL.mean():.6f} {chi:.6f}\n")
PY
done

echo "== analyze: deconfinement beta_c (the melting) =="
python3 - "$OUT/polyakov.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
beta, _, chi = np.loadtxt(sys.argv[1], unpack=True)
out = lat.deconfinement_beta_c(beta, chi)
print(f"  beta_c = {out['beta_c']:.4f}   (chi_L peak = {out['chi_peak']:.3f})")
print(f"  -> T_c = 1/(Nt * a(beta_c)); set a from target 3 (gradient flow).")
PY
