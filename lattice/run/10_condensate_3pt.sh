#!/usr/bin/env bash
# Target: the GENUINE in-medium condensate bag (Eternal Dawn, Ch. 11 Generations) -- the connected
# nucleon scalar 3-point <N|qbar q(r)|N> (measure_condensate_3pt), the definitive s_T behind the lever,
# refining measure_bag_profile's one-body proxy (run/09). For each config we solve the forward
# propagator, build the sequential source (the derivative of the validated nucleon contraction), solve
# the sequential propagator, and emit the scalar 3-point density binned by radius.
#
# VALIDATE FIRST. The 3-point needs the sequential source to be right. The program emits a SELF-CHECK:
# the source reconstruction MUST equal C_N(t_snk) (C_N is linear in the spectator propagator). The
# analysis checks this BEFORE reporting anything -- so the first job of this run is to confirm the gate
# passes at a clean heavy mass (default -0.5). Once it does, rerun at the light masses (-0.85, -0.95)
# where the bag is in/near the window to read the genuine condensate s_T. Cost: ~2x run/09 (a second
# 12-column solve for the sequential propagator).
set -euo pipefail
source "$(dirname "$0")/_lib.sh"
HERE="$(cd "$(dirname "$0")" && pwd)"
MEAS="$HERE/../src/measure_condensate_3pt"
require_exe "$MEAS"

OUT="${OUT:-out/dyn_L12x24_m-0.5}"
L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
MASS="${MASS:--0.5}"                    # start at a clean heavy mass to VALIDATE the self-check
# nucleon sink time slice (insertion scanned 0<tau<SINKT). MUST sit in the 2-pt plateau where C_N>0:
# a heavy nucleon decays fast, so T/2 can be PAST the plateau where the backward (wrong-parity) state
# flips C_N negative and contaminates the profile. ~3T/8 lands in the plateau here; the analysis warns
# if C_N(t_snk)<=0. Lower SINKT (heavier mass) / raise it (lighter) to track the 2-pt plateau.
SINKT="${SINKT:-$((3 * T / 8))}"
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-2}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
TAU_LO="${TAU_LO:-}" ; TAU_HI="${TAU_HI:-}"      # tau plateau (default: midway between source and sink)
R0A="${R0A:-3.166}"
# source SMEARING: widens the ground-state overlap so the 3-pt plateaus at small t_snk (a point
# source needs ~8 slices; smeared, far fewer). 0 = point (excited-state contaminated -- biases s_T
# LOW). radius ~ sqrt(2*iters*alpha); 20 x 0.25 ~ 3 lattice units, a nucleon-scale source.
SMEAR_ITERS="${SMEAR_ITERS:-20}" ; SMEAR_ALPHA="${SMEAR_ALPHA:-0.25}"

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
echo "condensate 3pt on ${#sel[@]} configs, mass=$MASS, sink-time=$SINKT, smear=${SMEAR_ITERS}x${SMEAR_ALPHA}, $NPAR in parallel"

start_mps
for cfg in "${sel[@]}"; do
  printf '%q --grid %q --config %q --mass %q --sink-time %q --cg-tol %q --smear-iters %q --smear-alpha %q --accelerator-threads 8 | grep -E "^(CHK|C2|SR|P3)" > %q\n' \
    "$MEAS" "$GRIDSPEC" "$cfg" "$MASS" "$SINKT" "$CGTOL" "$SMEAR_ITERS" "$SMEAR_ALPHA" "$cfg.c3pt"
done | run_pool "$NPAR"
stop_mps

ci=0; : > "$OUT/c3pt_raw.dat"
for cfg in "${sel[@]}"; do
  ci=$((ci+1))
  awk -v c=$ci '
    /^CHK/{print c, "CHK", $3, $4}
    /^C2/ {print c, "C2",  $2, $3}
    /^SR/ {print c, "SR",  $2, $3}
    /^P3/ {print c, "P3",  $2, $3, $4, $5}' "$cfg.c3pt" >> "$OUT/c3pt_raw.dat"
done
echo "measured ${#sel[@]} configs -> $OUT/c3pt_raw.dat"

echo "== analyze: the self-check, the scalar charge, and the condensate bag =="
python3 - "$OUT/c3pt_raw.dat" "$T" "$SINKT" "$R0A" "${TAU_LO:-_}" "${TAU_HI:-_}" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
rows = [l.split() for l in open(sys.argv[1]) if l.split()]
T, t_snk, r0a = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
tw = None
if sys.argv[5] != "_" and sys.argv[6] != "_":
    tw = (int(sys.argv[5]), int(sys.argv[6]))
chk = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if r[1] == "CHK"])
c2  = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if r[1] == "C2"])
sr  = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if r[1] == "SR"])
p3  = np.array([[float(x) for x in (r[0], r[2], r[3], r[4], r[5])] for r in rows if r[1] == "P3"])
res = lat.condensate_3pt(p3, c2, sr, chk, T=T, t_snk=t_snk, tau_window=tw, r0_over_a=r0a)

print(f"  SELF-CHECK: recon {res['recon']:.6g}  vs  C_N(t_snk) {res['cn_snk']:.6g}  "
      f"(resid {res['self_check_resid']:.2%})  -> {'PASS' if res['self_check_ok'] else 'FAIL'}")
print(f"  SINK     : C_N(t_snk={t_snk}) = {res['cn_snk']:.3e}  (2-pt turns non-positive at t={res['node_t']})"
      f"  -> {'OK' if res['sink_ok'] else 'TOO FAR (reduce SINKT)'}")
if res["self_check_ok"] and res["sink_ok"]:
    print(f"  scalar charge g_S = {res['g_S']:.3f}  (sum-rule plateau, tau in {res['tau_window']})")
    print(f"  condensate bag rho3(r):")
    for r, rr in zip(res["r"], res["rho3"]):
        bar = "#" * int(50 * abs(rr) / (max(abs(res['rho3'])) + 1e-300))
        print(f"    r={r:5.2f}  G3={rr:+10.4e}  {bar}")
    print(f"  half-density radius R0 = {res['R0']:.2f} a   s_T = R0/r0 = "
          f"{res['s_T']:.3f} +/- {res['s_T_err']:.3f}")
print(f"  VERDICT: {res['verdict']}")
PY
