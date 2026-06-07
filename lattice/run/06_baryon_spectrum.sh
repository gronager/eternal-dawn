#!/usr/bin/env bash
# Target L4 (pilot): the torsiton ground state -- quenched valence baryon spectroscopy on the
# pure-gauge ensemble. Reuses the L1 configs (no new HMC): solve a Wilson-quark propagator on each
# thermalised, decorrelated config, contract the PION (mass calibration / chiral probe) and the
# NUCLEON (the torsiton), and read the baryon mass off the effective-mass plateau.
#
# "Rediscover the bound ground state on a lattice": a clean plateau at m_N > 0 IS the torsiton,
# found non-perturbatively. Cheap pilot; dynamical fermions near the chiral limit (the real
# spectrum, and any excited rungs = candidate further generations) come after.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, MPS, run_pool
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/puregauge}"                        # SAME ensemble as 01_puregauge_potential.sh
MEAS="$HERE/../src/measure_baryon"
require_exe "$MEAS"

L="${L:-16}" ; T="${T:-32}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:-0.1}"        # Wilson bare quark mass. Heavier = cleaner signal; lighter = more chiral.
                          # Scan a few (e.g. 0.2, 0.1, 0.05) to extrapolate toward the chiral limit.
THERM="${THERM:-150}"     # thermalisation cut (trajectory index)
STRIDE="${STRIDE:-20}"    # decorrelation stride
NPAR="${NPAR:-8}"         # parallel solves on the GPU (drop toward 1-2 once one lattice fills it)
CGTOL="${CGTOL:-1e-8}"

# --- select thermalised, decorrelated configs (same logic as 01) ---------------------------
shopt -s nullglob
sel=()
for cfg in "$OUT"/ckpoint_lat.* "$OUT"/stream*/ckpoint_lat.*; do
  [[ "$cfg" == *.gz ]] && continue
  n="${cfg##*.}"; [[ "$n" =~ ^[0-9]+$ ]] || continue
  (( n >= THERM )) || continue
  (( n % STRIDE == 0 )) || continue
  sel+=("$cfg")
done
[ "${#sel[@]}" -gt 0 ] || { echo "no thermalised configs in $OUT (run 00/01 first)"; exit 1; }
echo "measuring baryon+pion on ${#sel[@]} configs, mass=$MASS, $NPAR in parallel"

# --- solve + contract on each config (parallel on the GPU) ---------------------------------
start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --cg-tol %q --accelerator-threads 8 | grep -E "^[0-9]" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$CGTOL" "$cfg.bar"
done | run_pool "$NPAR"
stop_mps

# combine, tagging each row with a config index (for the jackknife)
ci=0; : > "$OUT/baryon_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '/^[0-9]/{print c, $0}' "$cfg.bar" >> "$OUT/baryon_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/baryon_raw.dat"

# --- analyse: effective masses + plateau fits (jackknife over configs) ---------------------
echo "== analyze: pion and torsiton (nucleon) masses =="
python3 - "$OUT/baryon_raw.dat" "$T" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split() and l.split()[0].lstrip("-").isdigit()]
raw = np.array([[float(x) for x in r] for r in rows if len(r) == 4])   # cfg t C_pi C_N
if raw.size == 0:
    sys.exit("no 'cfg t C_pi C_N' rows in baryon_raw.dat")
res = lat.baryon_spectrum(raw, T=int(sys.argv[2]))
for name in ("pion", "nucleon"):
    r = res[name]
    print(f"\n  {name}: effective mass m_eff(t)=ln[C(t)/C(t+1)]")
    print("    " + "  ".join(f"{t}:{m:.3f}" for t, m in zip(r['t'], r['meff']) if np.isfinite(m)))
    if np.isfinite(r['mass']):
        print(f"    plateau fit  m = {r['mass']:.4f} +/- {r['mass_err']:.4f}  "
              f"(window t=[{r['tmin']},{r['tmax']}], n_cfg={r['n_cfg']})")
    else:
        print("    no clean plateau -- adjust the window / raise statistics / heavier mass")
if np.isfinite(res['nucleon']['mass']) and np.isfinite(res['pion']['mass']):
    print(f"\n  m_N / m_pi = {res['nucleon']['mass']/res['pion']['mass']:.3f}   "
          f"(the torsiton is bound at m_N>0; ratio drops toward the chiral limit)")
PY
