#!/usr/bin/env bash
# Target: the s_T lever (Eternal Dawn, Ch. 11 Generations). Measure the SHAPE of the torsiton's
# mass-giving bag -- the one number the fermion-mass hierarchy reduces to. For each thermalised
# config we solve a point-source propagator and bin its gauge-invariant scalar density
# rho(r)=Tr[S(x;0)^dag S(x;0)] by spatial radius on each time slice (measure_bag_profile). The
# config-averaged rho(r) in the plateau is the bag: its half-density radius R0 and wall thickness a
# fix the sharpness s_T = R0/r0, which the analysis crosses against the lepton lever
# (cartasis_sims.fermion_masses). The lever reproduces 3477 only for s_T in a productive window
# [0.43,0.70] r0; the verdict flags whether the measured bag lands there (heavy pions make it small).
#
# Runs on the dynamical ensemble (set OUT to the run/07 output). Cheap: one propagator solve per
# config, no sink contraction -- packs many in parallel under MPS like run/06.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"
HERE="$(cd "$(dirname "$0")" && pwd)"
MEAS="$HERE/../src/measure_bag_profile"
require_exe "$MEAS"

OUT="${OUT:-out/dyn_L12x24_m-0.5}"
L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:--0.5}"                  # valence = sea (unitary point), as run/07
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-2}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
# plateau window for averaging rho(r,t) (pion plateau on this ensemble ~ t in [4,10]); analysis default
PLAT_LO="${PLAT_LO:-4}" ; PLAT_HI="${PLAT_HI:-10}"

shopt -s nullglob
sel=()
for cfg in "$OUT"/ckpoint_lat.* "$OUT"/stream*/ckpoint_lat.*; do
  [[ "$cfg" == *.gz ]] && continue
  n="${cfg##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
  (( n >= THERM )) || continue
  (( n % STRIDE == 0 )) || continue
  sel+=("$cfg")
done
[ "${#sel[@]}" -gt 0 ] || { echo "no configs in $OUT (n>=$THERM, stride $STRIDE)"; exit 1; }
echo "bag profile on ${#sel[@]} configs, mass=$MASS, $NPAR in parallel"

start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --cg-tol %q --accelerator-threads 8 | grep -E "^(PROF|PION)" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$CGTOL" "$cfg.bag"
done | run_pool "$NPAR"
stop_mps

ci=0; : > "$OUT/bag_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '/^PROF/{print c, $2, $3, $4, $5} /^PION/{print c, "PION", $2, $3}' "$cfg.bag" >> "$OUT/bag_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/bag_raw.dat"

echo "== analyze: the bag shape s_T and the lepton lever =="
python3 - "$OUT/bag_raw.dat" "$T" "$PLAT_LO" "$PLAT_HI" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split()]
T, lo, hi = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
prof = np.array([[float(x) for x in r] for r in rows if r[1] != "PION" and len(r) == 5])  # cfg r2 t sum cnt
res = lat.bag_profile(prof, T=T, plateau=(lo, hi))
print(f"  config-averaged bag rho(r), plateau t in [{lo},{hi}], n_cfg={res['n_cfg']}")
for r, rr in zip(res["r"], res["rho"]):
    bar = "#" * int(60 * rr / max(res["rho"]))
    print(f"    r={r:5.2f}  rho={rr:10.4e}  {bar}")
print(f"  half-density radius R0 = {res['R0']:.2f} a   (Fermi wall thickness {res['wall']:.2f} a)")
print(f"  sharpness  s_T = R0/r0 = {res['s_T']:.3f} +/- {res['s_T_err']:.3f}   (r0/a = {res['r0_over_a']:.3f})")
print(f"  -> lever span at this s_T: {res['span']:.0f} +/- {res['span_err']:.0f}   "
      f"(observed 3477; trustworthy only for s_T in {res['productive_window']} r0)")
print(f"  VERDICT: {res['verdict']}")
PY
