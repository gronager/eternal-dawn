#!/usr/bin/env bash
# Target L3 (the control): a FIRST lattice estimate of the electroweak S parameter (Peskin-Takeuchi)
# on the dynamical N_f=2 SU(3) FUNDAMENTAL Wilson ensemble, via the isovector VECTOR minus AXIAL-VECTOR
# correlators (Pi_V - Pi_A). This tests the framework's load-bearing assumption (Eternal Dawn, Ch. "The
# Forces from the Field"): does the electroweak-breaking sector clear S < 0.1, or sit in the QCD-like
# graveyard at S ~ 0.25?
#
# For each config, solve a point-source Wilson-quark propagator S and build the connected isovector
# correlators C_V(t) (gamma_i, rho) and C_A(t) (gamma_i*gamma5, a1), summed over the spatial
# polarisations i=1,2,3 (measure_s_parameter). The analysis cosh-fits C_V and C_A for M_V (rho), M_A
# (a1), and the decay constants F_V, F_A, then forms the lowest-resonance (VMD / Weinberg-sum-rule-
# saturated) proxy S = 4*pi*(F_V^2/M_V^2 - F_A^2/M_A^2), and prints M_A/M_V (the walking diagnostic)
# and the integrated V-A (the chiral order parameter).
#
# HONEST: the FUNDAMENTAL N_f=2 rep is QCD-like (NOT near-conformal), so this proxy is EXPECTED to land
# near the QCD value S ~ 0.25 -- that CHARACTERISES the sector and TESTS the assumption, it does NOT
# validate S < 0.1 (which needs the near-conformal / propagating-torsion regime, target L5). The
# currents are UNRENORMALISED (Z_V, Z_A owed) at a single heavy sea mass. Do not oversell.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, MPS, run_pool, PYTHONPATH
HERE="$(cd "$(dirname "$0")" && pwd)"
MEAS="$HERE/../src/measure_s_parameter"
require_exe "$MEAS"

OUT="${OUT:-out/dyn_L12x24_m-0.5}"                 # the dynamical N_f=2 fundamental ensemble
L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:--0.5}"        # Wilson bare INPUT mass of the ensemble (the sea/valence mass)
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-2}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
TMIN="${TMIN:-}" ; TMAX="${TMAX:-}"                # cosh fit window (default: T//8 .. T//2-1)

# --- select thermalised, decorrelated configs (same logic as 06/10/11) ----------------------
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
echo "S parameter (V-A) on ${#sel[@]} configs, mass=$MASS, $NPAR in parallel"

# --- solve + contract the V and A correlators on each config (parallel on the GPU) -----------
start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --cg-tol %q --accelerator-threads 8 | grep -E "^(V|A) " > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$CGTOL" "$cfg.va"
done | run_pool "$NPAR"
stop_mps

# combine, tagging each row with a config index (for the jackknife): cfg CH t C(t)
ci=0; : > "$OUT/sparam_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '/^V /{print c, "V", $2, $3} /^A /{print c, "A", $2, $3}' "$cfg.va" >> "$OUT/sparam_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/sparam_raw.dat"

# --- analyse: cosh fits M_V, M_A, F_V, F_A; the S-proxy; M_A/M_V; the integrated V-A verdict ---
echo "== analyze: S = 4*pi*(F_V^2/M_V^2 - F_A^2/M_A^2), the QCD-like control =="
python3 - "$OUT/sparam_raw.dat" "$T" "${TMIN:-_}" "${TMAX:-_}" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split()]
T = int(sys.argv[2])
tmin = None if sys.argv[3] == "_" else int(sys.argv[3])
tmax = None if sys.argv[4] == "_" else int(sys.argv[4])
c_v = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if len(r) == 4 and r[1] == "V"])
c_a = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if len(r) == 4 and r[1] == "A"])
if c_v.size == 0 or c_a.size == 0:
    sys.exit("no 'cfg V/A t C(t)' rows in sparam_raw.dat")
res = lat.s_parameter_proxy(c_v, c_a, T=T, tmin=tmin, tmax=tmax)

print(f"  fit window t=[{res['tmin']},{res['tmax']}], n_cfg={res['n_cfg']}, npol={res['npol']}")
print(f"  M_V*a (rho) = {res['M_V']:.4f} +/- {res['M_V_err']:.4f}")
print(f"  M_A*a (a1)  = {res['M_A']:.4f} +/- {res['M_A_err']:.4f}")
print(f"  M_A/M_V     = {res['ratio']:.3f} +/- {res['ratio_err']:.3f}   (->1 conformal, ~1.6 QCD-like)")
print(f"  F_V*a       = {res['F_V']:.4e} +/- {res['F_V_err']:.2e}   (UNRENORMALISED, Z_V owed)")
print(f"  F_A*a       = {res['F_A']:.4e} +/- {res['F_A_err']:.2e}   (UNRENORMALISED, Z_A owed)")
print(f"  sum_t (C_V - C_A) = {res['va_int']:.4e} +/- {res['va_int_err']:.2e}   (chiral order parameter)")
print(f"  S-proxy = {res['S']:.4f} +/- {res['S_err']:.4f}   (vs QCD ~{res['S_qcd_ref']:.2f}, target S<{res['S_target']:.2f})")
print(f"  VERDICT: {res['verdict']}")
PY
