#!/usr/bin/env bash
# RUNG 3 -- THE GATE: gamma_m of SU(3) N_f=2 sextet from the Dirac mode number, on a DYNAMICAL
# ensemble (sextet SEA quarks). Two stages: (a) generate the ensemble with our generate_sextet HMC
# (Grid two-index HMC, CLI beta/mass); (b) measure the sextet mode number on the configs with the
# validated KPM chain and fit gamma_m. The verdict: nu(M) ~ M^{4/(1+gamma_m)} in the scaling window
# -> gamma_m > ~0.15 favours a walking sector (Part III), QCD-like/conformal retreats to Part II.
#
# (beta, mass) ARE THE PHYSICS. Bracket the critical mass first with run/04 (quenched), then tune
# here: mass toward critical, beta in the scaling window. Expect days-weeks on a single GH200.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, PYTHONPATH, MPS, run_pool
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/sextet}"; mkdir -p "$OUT"
GEN="$HERE/../src/generate_sextet"
MEAS="$HERE/../src/measure_modenumber"
require_exe "$GEN"
require_exe "$MEAS"

L="${L:-16}" ; T="${T:-32}" ; GRIDSPEC="$L.$L.$L.$T"
BETA="${BETA:-5.4}" ; MASS="${MASS:--0.95}"         # TUNE these (see run/04 for the m_c bracket)
NTRAJ="${NTRAJ:-2000}" ; THERM="${THERM:-400}" ; STRIDE="${STRIDE:-20}"
MDSTEPS="${MDSTEPS:-20}" ; SAVE_INTERVAL="${SAVE_INTERVAL:-5}"
NMOM="${NMOM:-800}" ; NNOISE="${NNOISE:-12}" ; XMAX="${XMAX:-75}"
MLO="${MLO:-0.04}" ; MHI="${MHI:-0.30}"             # mode-number scaling window (drop lattice-artefact ends)

# --- (a) generate the dynamical sextet ensemble (skip if configs already exist) -------------
shopt -s nullglob
if ls "$OUT"/ckpoint_lat.* >/dev/null 2>&1; then
  echo "using existing sextet configs in $OUT"
else
  echo "== (a) dynamical HMC: SU(3) Nf=2 sextet  beta=$BETA mass=$MASS  $GRIDSPEC =="
  ( cd "$OUT" && "$GEN" --grid "$GRIDSPEC" --beta "$BETA" --mass "$MASS" --seed 1001 \
        --mdsteps "$MDSTEPS" --save-interval "$SAVE_INTERVAL" \
        --StartingType HotStart --Trajectories "$NTRAJ" --accelerator-threads 8 \
        2>&1 | tee hmc.log )
fi

# --- (b) sextet mode number on thermalised, decorrelated configs (skip .gz/.wl) -------------
sel=()
for cfg in "$OUT"/ckpoint_lat.*; do
  case "$cfg" in *.gz|*.wl) continue;; esac
  n="${cfg##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
  (( n >= THERM )) || continue
  (( n % STRIDE == 0 )) || continue
  sel+=("$cfg")
done
[ "${#sel[@]}" -gt 0 ] || { echo "no thermalised configs (n>=$THERM, stride $STRIDE) -- raise NTRAJ"; exit 1; }
echo "== (b) sextet mode number on ${#sel[@]} configs =="

mfiles=()
for cfg in "${sel[@]}"; do
  mf="$OUT/mom.$(basename "$cfg").txt"
  if "$MEAS" --grid "$GRIDSPEC" --config "$cfg" --rep sextet --mass "$MASS" \
             --xmax "$XMAX" --Nmom "$NMOM" --Nnoise "$NNOISE" > "$mf" 2>"$mf.log"; then
    mfiles+=("$mf"); echo "  measured $(basename "$cfg")"
  else
    echo "  measure failed ($(basename "$cfg")) -- see $mf.log"
  fi
done
[ "${#mfiles[@]}" -gt 0 ] || { echo "no moment files produced"; exit 1; }

echo "== gamma_m (THE GATE) =="
"${PY:-python3}" - "$MLO" "$MHI" "${mfiles[@]}" <<'PY' | tee "$OUT/gamma_m.dat"
import sys, numpy as np
from cartasis_sims import lattice as lat
mlo, mhi = float(sys.argv[1]), float(sys.argv[2]); files = sys.argv[3:]
M = np.linspace(mlo * 0.7, mhi * 1.3, 50); nus = []
for f in files:
    xm = None; mu = {}
    for ln in open(f):
        t = ln.split()
        if len(t) >= 2 and t[0] == "XMAX": xm = float(t[1])
        elif len(t) >= 3 and t[0] == "MOMENT": mu[int(t[1])] = float(t[2])
    if xm and mu:
        nus.append(lat.mode_number_from_chebyshev_moments(
            np.array([mu[k] for k in range(max(mu) + 1)]), M, xm))
nu = np.mean(nus, axis=0)
o = lat.anomalous_dimension_from_mode_number(M, nu, window=(mlo, mhi))
g = o["gamma_m"]
print(f"  configs = {len(nus)}   window M in [{mlo}, {mhi}]")
print(f"  gamma_m = {g:+.3f}   (slope {o['slope']:.3f})")
print("  VERDICT:", "WALKS -> Part III favourably placed (chase S)" if g > 0.15
      else "QCD-like / conformal-window -> Part III retreats to Part II")
PY
echo "-> $OUT/gamma_m.dat"
