#!/usr/bin/env bash
# Target L4 (excited spectrum): the EXCITED torsiton tower by the variational method (GEVP) -- the
# direct lattice test of the three-generation thesis. The framework predicts exactly three radial
# rungs (the generations) and no fourth. A single operator (run/06) sees only the ground state; here
# we build an N-operator basis (the nucleon at N Gaussian smearing radii), form the correlator
# matrix C_ij(t), and solve C(t) v = lambda(t) C(t0) v -- whose N eigenvalues resolve N states.
#
# Reads VERDICT: count the clean levels. If a third nucleon state plateaus and a fourth is noise,
# that supports three rungs; if the tower keeps going, it does not. Works on quenched OR dynamical
# configs (set OUT). Excited-state spectroscopy is hard -- expect the ground state clean, the first
# excited fair, the higher ones noisy on a small/coarse pilot; the signal sharpens with statistics,
# more operators, and (for the radial structure) a bigger box.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/dyn_L12x24_m-0.5}"
MEAS="$HERE/../src/measure_baryon_gevp"
require_exe "$MEAS"

L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:--0.5}"
ITERS="${ITERS:-0,20,50}"      # Wuppertal smearing-step counts; one operator each; 0 = point. radius ~ sqrt(2*iters*alpha)
SALPHA="${SALPHA:-0.25}"       # smearing weight (normalised form -- stable for any alpha)
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-10}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
NOP=$(awk -F, '{print NF}' <<<"$ITERS")

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
echo "GEVP on ${#sel[@]} configs, $NOP operators (smear-iters $ITERS), mass=$MASS, $NPAR in parallel"

start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --smear-iters %q --smear-alpha %q --cg-tol %q --accelerator-threads 8 | grep -E "^[0-9]" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$ITERS" "$SALPHA" "$CGTOL" "$cfg.gevp"
done | run_pool "$NPAR"
stop_mps

ci=0; : > "$OUT/gevp_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '/^[0-9]/{print c, $0}' "$cfg.gevp" >> "$OUT/gevp_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/gevp_raw.dat"

echo "== analyze: the excited torsiton tower (GEVP) =="
python3 - "$OUT/gevp_raw.dat" "$NOP" "$T" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split() and l.split()[0].lstrip("-").isdigit()]
raw = np.array([[float(x) for x in r] for r in rows if len(r) == 5])   # cfg i j t C
N, T = int(sys.argv[2]), int(sys.argv[3])
if raw.size == 0:
    sys.exit("no 'cfg i j t C' rows in gevp_raw.dat")
res = lat.gevp_spectrum(raw, N=N, T=T)
print(f"  {N}-operator GEVP (t0={res['t0']}, plateau t=[{res['tmin']},{res['tmax']}], n_cfg={res['n_cfg']})")
names = ["ground (gen I)", "1st excited (gen II)", "2nd excited (gen III)", "3rd excited (gen IV?)"]
for n in range(N):
    m, e = res["masses"][n], res["mass_err"][n]
    tag = names[n] if n < len(names) else f"level {n}"
    clean = "clean" if (np.isfinite(e) and e < 0.15 * abs(m)) else "NOISY (basis running out?)"
    print(f"    {tag:24s}: m = {m:.4f} +/- {e:.4f}   [{clean}]")
fin = np.isfinite(res["masses"])
if fin.sum() >= 2:
    g = res["masses"][0]
    print("  gaps above ground:", "  ".join(f"{(res['masses'][n]-g):.3f}" for n in range(1, N) if fin[n]))
print("  VERDICT: count the CLEAN levels -- three clean rungs + a noisy fourth supports the thesis;")
print("           a clean fourth (or no clean second) does not. Sharpen with stats/operators/box.")
PY
