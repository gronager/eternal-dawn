#!/usr/bin/env bash
# Quenched mode-number scan: gamma_m(mass) for a representation on the pure-gauge ensemble.
# This is a VALIDATION / PREP rung, not the physics gamma_m. It sweeps the Wilson mass and watches
# the mode-number slope: far from criticality the spectrum is gapped (slope steep / nu empty), and
# as the mass approaches the critical point the low modes reach toward zero and the slope settles
# near 4/(1+gamma_m). Two jobs: (a) confirm the sextet operator behaves PHYSICALLY vs mass on real
# links (the link map is Grid's update_representation; 2a only tested identity links), and (b)
# bracket the sextet critical mass as a prior for the dynamical run. The walks-vs-conformal verdict
# needs the DYNAMICAL sea (05), not this quenched valence scan.
#
# Env: REP, MASSES, MAXCFG, MLO/MHI, NMOM/NNOISE/XMAX, IN/OUT. The analysis uses numpy + cartasis_sims;
# if the box python3 lacks numpy, point PY at one that has it, e.g. PY=/tmp/ed-venv/bin/python ...
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, PYTHONPATH (cartasis_sims)
HERE="$(cd "$(dirname "$0")" && pwd)"
IN="${IN:-out/puregauge}"                           # where 00/01 put the pure-gauge configs
OUT="${OUT:-out/quenched_modenumber}"; mkdir -p "$OUT"
MEAS="$HERE/../src/measure_modenumber"
require_exe "$MEAS"

REP="${REP:-sextet}"                                # fundamental | sextet
MASSES="${MASSES:--0.6 -0.4 -0.2 0.0 0.2 0.4 0.6}"  # Wilson mass sweep
NMOM="${NMOM:-800}" ; NNOISE="${NNOISE:-12}" ; XMAX="${XMAX:-75}"   # 75 bounds Wilson D^dagD at |m|<=8
THERM="${THERM:-150}" ; STRIDE="${STRIDE:-20}" ; MAXCFG="${MAXCFG:-2}"
MLO="${MLO:-0.20}" ; MHI="${MHI:-0.45}"

# --- collect thermalised, decorrelated configs (skip .gz and .wl side-products), cap to MAXCFG ---
shopt -s nullglob
sel=()
for cfg in "$IN"/ckpoint_lat.* "$IN"/stream*/ckpoint_lat.*; do
  case "$cfg" in *.gz|*.wl) continue;; esac
  n="${cfg##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
  (( n >= THERM )) || continue
  (( n % STRIDE == 0 )) || continue
  sel+=("$cfg")
done
[ "${#sel[@]}" -gt 0 ] || { echo "no configs in $IN (n>=$THERM, stride $STRIDE)"; exit 1; }
(( ${#sel[@]} > MAXCFG )) && sel=("${sel[@]: -MAXCFG}")    # newest MAXCFG

# --- grid from the NERSC header (authoritative; no guessing) ---
GRIDSPEC="$(head -c 4000 "${sel[0]}" | awk -F'[[:space:]=]+' \
  '/DIMENSION_1/{a=$2}/DIMENSION_2/{b=$2}/DIMENSION_3/{c=$2}/DIMENSION_4/{d=$2}END{print a"."b"."c"."d}')"
[ -n "$GRIDSPEC" ] || { echo "could not read DIMENSION_* from ${sel[0]}"; exit 1; }
echo "rep=$REP  grid=$GRIDSPEC  configs=${#sel[@]}  window=[$MLO,$MHI]  masses: $MASSES"

printf '# %-8s %-12s %-8s\n' mass gamma_m slope | tee "$OUT/gamma_vs_mass.dat"
for m in $MASSES; do
  mfiles=()
  for cfg in "${sel[@]}"; do
    mf="$OUT/mom.$REP.m$m.$(basename "$cfg").txt"
    if "$MEAS" --grid "$GRIDSPEC" --config "$cfg" --rep "$REP" --mass "$m" \
               --xmax "$XMAX" --Nmom "$NMOM" --Nnoise "$NNOISE" > "$mf" 2>"$mf.log"; then
      mfiles+=("$mf")
    else
      echo "  measure failed (mass $m, $(basename "$cfg")) -- see $mf.log"
    fi
  done
  [ "${#mfiles[@]}" -gt 0 ] || continue
  "${PY:-python3}" - "$MLO" "$MHI" "$m" "${mfiles[@]}" <<'PY' | tee -a "$OUT/gamma_vs_mass.dat"
import sys, numpy as np
from cartasis_sims import lattice as lat
mlo, mhi, mass = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
files = sys.argv[4:]
M = np.linspace(mlo * 0.8, mhi * 1.2, 40)
nus = []
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
print(f"  {mass:<8} {o['gamma_m']:<+12.4f} {o['slope']:<8.3f}")
PY
done
echo "-> $OUT/gamma_vs_mass.dat   (look for slope -> ~4 / gamma_m -> small as the spectrum nears criticality)"
