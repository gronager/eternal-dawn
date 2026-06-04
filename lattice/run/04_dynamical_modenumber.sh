#!/usr/bin/env bash
# Target 4 (THE GATE): the mass anomalous dimension gamma_m of the candidate walking sector,
# from the Dirac MODE NUMBER nu(M) ~ M^{4/(1+gamma_m)} (Giusti-Luscher, stochastic estimator).
# Candidate: SU(3), N_f=2 sextet (2-index symmetric). Cost: days-weeks, single GH200, no comms.
#
# Two stages: (a) generate a dynamical-fermion ensemble (HMC with the sextet Dirac operator);
# (b) on each config, stochastically count Dirac eigenvalues below a ladder of thresholds M.
set -euo pipefail
GRID="${GRID:-$HOME/ed-lattice/src/Grid/build}"
OUT="${OUT:-out/gammam}"; mkdir -p "$OUT"

L=24 ; T=48
BETA=5.4 ; MASS=-0.02      # bare coupling + (near-critical) fermion mass: TUNE these first
NCFG=300 ; NSEP=10
MGRID="0.02 0.04 0.06 0.08 0.10 0.14 0.18 0.24 0.30 0.40"   # thresholds for the mode number

echo "== (a) dynamical HMC: SU(3) N_f=2 sextet =="
# NOTE: requires a Grid build with the higher-rep (2-index-symmetric) fermion action.
"$GRID/tests/hmc/Test_hmc_WilsonFermionGauge" \
    --grid $L.$L.$L.$T --beta $BETA --mass $MASS --representation sextet --Nf 2 \
    --Trajectories $((NCFG*NSEP)) --Thermalizations 200 \
    --Checkpoint "$OUT/cfg" --accelerator-threads 8  2>&1 | tee "$OUT/hmc.log"

echo "== (b) Dirac mode number nu(M) (stochastic; no eigenvectors stored) =="
: > "$OUT/modenumber.dat"   # columns: M  nu  nu_err
for M in $MGRID; do
  # stochastic trace of the projector onto eigenmodes below M, averaged over configs+sources
  "$GRID/tests/solver/Test_ModeNumber" --grid $L.$L.$L.$T --representation sextet \
      --mass $MASS --threshold "$M" --nstoch 12 --configs "$OUT"/cfg.* \
      >> "$OUT/modenumber.dat"
done

echo "== analyze: gamma_m (THE GATE) =="
python3 - "$OUT/modenumber.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
M, nu = np.loadtxt(sys.argv[1], usecols=(0,1), unpack=True)
# fit in the scaling window (drop the lowest/highest M where lattice artefacts bite)
out = lat.anomalous_dimension_from_mode_number(M, nu, window=(0.04, 0.30))
g = out['gamma_m']
print(f"  gamma_m = {g:.3f}   (slope {out['slope']:.3f})")
print("  VERDICT:", "WALKS -> Part III favourably placed (chase S)" if g > 0.15
      else "QCD-like / conformal-window -> Part III retreats to Part II")
PY
