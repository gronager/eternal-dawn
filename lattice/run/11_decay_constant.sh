#!/usr/bin/env bash
# Target L2: the substrate scale v (= f_pi) -- the missing input for the framework's SHARP prediction
#
#     sigma = 2 * pi * v^2     (ONE substrate sets BOTH the mass scale and the confinement tension).
#
# We already have sqrt(sigma)*a = 0.385 (static potential, run/01) and m_pi*a = 1.269 (spectrum,
# run/06). This run supplies f_pi*a on the EXISTING dynamical N_f=2 SU(3) FUNDAMENTAL Wilson ensemble:
# for each config, solve a point-source Wilson-quark propagator S and build the zero-momentum
# pseudoscalar correlator C_PP(t) = sum_x Tr[S^dag S] (measure_decay_constant). The analysis fits the
# two-sided cosh for m_pi*a and the pseudoscalar overlap G_PP, then takes f_pi*a via the PCAC relation
# f_pi*a = (2 m_q / m_pi^2) sqrt(G_PP), and prints the sigma = 2*pi*v^2 verdict (v = f_pi, the
# framework's NJL identification).
#
# HONEST: f_pi*a here is UNRENORMALISED -- Z_A, Z_P and the critical mass m_crit are all OWED (each an
# O(1) factor). m_q = MASS - MCRIT (set MCRIT once it is measured; default 0 reports the bare-input
# m_q). Convention: f_pi ~ 92 MeV (F_pi); x sqrt(2) for the 130 MeV convention.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, MPS, run_pool, PYTHONPATH
HERE="$(cd "$(dirname "$0")" && pwd)"
MEAS="$HERE/../src/measure_decay_constant"
require_exe "$MEAS"

OUT="${OUT:-out/dyn_L12x24_m-0.5}"                 # the dynamical N_f=2 fundamental ensemble
L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:--0.5}"        # Wilson bare INPUT mass of the ensemble (the sea/valence mass)
MCRIT="${MCRIT:-0}"         # critical mass (OWED): m_q = MASS - MCRIT. Default 0 -> report bare m_q.
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-2}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
TMIN="${TMIN:-}" ; TMAX="${TMAX:-}"                # cosh fit window (default: T//8 .. T//2-1)
SQSIGMA="${SQSIGMA:-0.385}"                        # sqrt(sigma)*a from the static potential (run/01)

# --- select thermalised, decorrelated configs (same logic as 06/10) -------------------------
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
echo "decay constant on ${#sel[@]} configs, mass=$MASS, mcrit=$MCRIT, $NPAR in parallel"

# --- solve + contract the PP correlator on each config (parallel on the GPU) -----------------
start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --cg-tol %q --accelerator-threads 8 | grep -E "^PP" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$CGTOL" "$cfg.pp"
done | run_pool "$NPAR"
stop_mps

# combine, tagging each row with a config index (for the jackknife): cfg t C_PP(t)
ci=0; : > "$OUT/decay_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '/^PP/{print c, $2, $3}' "$cfg.pp" >> "$OUT/decay_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/decay_raw.dat"

# --- analyse: cosh fit for m_pi*a and G_PP, f_pi*a (PCAC), and the sigma = 2*pi*v^2 verdict ----
echo "== analyze: f_pi*a (= v*a) and the sigma = 2*pi*v^2 test =="
MQ=$(python3 -c "print($MASS - $MCRIT)")
python3 - "$OUT/decay_raw.dat" "$T" "$MQ" "$SQSIGMA" "${TMIN:-_}" "${TMAX:-_}" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split() and l.split()[0].lstrip("-").isdigit()]
raw = np.array([[float(x) for x in r] for r in rows if len(r) == 3])     # cfg t C_PP
if raw.size == 0:
    sys.exit("no 'cfg t C_PP' rows in decay_raw.dat")
T, m_q, sqsig = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
tmin = None if sys.argv[5] == "_" else int(sys.argv[5])
tmax = None if sys.argv[6] == "_" else int(sys.argv[6])
res = lat.decay_constant(raw, m_q=m_q, T=T, tmin=tmin, tmax=tmax)

print(f"  PCAC route: f_pi*a = (2 m_q / m_pi^2) sqrt(G_PP), m_q={m_q:.4f} (bare-input - m_crit; m_crit OWED)")
print(f"  fit window t=[{res['tmin']},{res['tmax']}], n_cfg={res['n_cfg']}, {res['convention']}")
print(f"  m_pi*a  = {res['m_pi']:.4f} +/- {res['m_pi_err']:.4f}")
print(f"  G_PP    = {res['G_PP']:.4e} +/- {res['G_PP_err']:.2e}")
print(f"  f_pi*a  = {res['f_pi']:.4f} +/- {res['f_pi_err']:.4f}   (= v*a; UNRENORMALISED, Z_A/Z_P owed)")
chk = lat.sigma_2piv2_check(sqsig, res['f_pi'],
                            sqrt_sigma_a_err=None, f_pi_a_err=res['f_pi_err'])
print(f"  sqrt(sigma)*a = {sqsig:.3f} (run/01)  ->  predicted f_pi*a = {chk['f_pi_pred_a']:.4f}")
print(f"  ratio sigma/(2*pi*f_pi^2) = {chk['ratio']:.3f}"
      + (f" +/- {chk['ratio_err']:.3f}" if np.isfinite(chk['ratio_err']) else ""))
print(f"  VERDICT: {chk['verdict']}")
PY
