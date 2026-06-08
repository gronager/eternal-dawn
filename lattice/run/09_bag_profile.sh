#!/usr/bin/env bash
# Target: the s_T lever (Eternal Dawn, Ch. 11 Generations). Measure the SHAPE of the torsiton's
# mass-giving bag -- the one number the fermion-mass hierarchy reduces to -- AND its chiral trend.
# For each thermalised config and each valence quark mass we solve a point-source propagator and bin
# its gauge-invariant scalar density rho(r)=Tr[S(x;0)^dag S(x;0)] by spatial radius on each time slice
# (measure_bag_profile). The config-averaged rho(r) in the plateau is the bag: its half-density radius
# R0 fixes s_T = R0/r0; the spatial sum is the pion correlator, so each mass also gives m_pi.
#
# THE REAL TEST is not one heavy point but the TREND: a heavy quark gives a SMALL bag; does s_T RISE
# into the lever's productive window [0.43,0.70] r0 as m_pi^2 -> 0 (where the consecutive ladder
# reproduces the observed lepton span ~3477)? We scan a few VALENCE masses (partially quenched, on the
# run/07 sea ensemble -- no new HMC, cheap) and extrapolate s_T vs m_pi^2 to the chiral point.
#
# Cheap: one propagator solve per (config, mass), no sink contraction -- packs in parallel under MPS.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"
HERE="$(cd "$(dirname "$0")" && pwd)"
MEAS="$HERE/../src/measure_bag_profile"
require_exe "$MEAS"

OUT="${OUT:-out/dyn_L12x24_m-0.5}"
L="${L:-12}" ; T="${T:-24}" ; GRIDSPEC="$L.$L.$L.$T"
# valence masses (partially quenched), heavy -> light toward m_crit ~ -1.15; stop BEFORE the wall
# (m <~ -1.0): below it the valence quark goes critical and exceptional configs delocalise the bag to
# the whole box (m_pi -> 0, R0 -> L/2). A single value reproduces the old single-point behaviour.
VALMASS="${VALMASS:--0.5 -0.7 -0.85 -0.95}"
THERM="${THERM:-40}" ; STRIDE="${STRIDE:-2}" ; NPAR="${NPAR:-4}" ; CGTOL="${CGTOL:-1e-8}"
PLAT_LO="${PLAT_LO:-4}" ; PLAT_HI="${PLAT_HI:-10}"      # rho(r,t) / pion plateau window
R0A="${R0A:-3.166}"                                     # r0/a on this ensemble (run/01 dynamical)
# quality guard: drop exceptional points (delocalised bag or near-massless pion) from the trend
MAXR0FRAC="${MAXR0FRAC:-0.30}"                          # drop if R0 > MAXR0FRAC*L (bag fills the box)
MIN_MPI="${MIN_MPI:-0.10}"                              # drop if m_pi a < MIN_MPI (exceptional zero mode)

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
echo "bag profile on ${#sel[@]} configs, valence masses: $VALMASS  ($NPAR in parallel)"

# --- measure at each valence mass; collect one raw file per mass -----------------------------------
py_args=()
for m in $VALMASS; do
  tag="$(echo "$m" | sed 's/-/m/; s/\./p/')"          # -0.85 -> m0p85
  raw="$OUT/bag_raw.$tag.dat"
  echo "-- valence mass $m  ->  $raw"
  start_mps
  for cfg in "${sel[@]}"; do
    printf '%q --grid %q --config %q --mass %q --cg-tol %q --accelerator-threads 8 | grep -E "^(PROF|PION)" > %q\n' \
      "$MEAS" "$GRIDSPEC" "$cfg" "$m" "$CGTOL" "$cfg.bag.$tag"
  done | run_pool "$NPAR"
  stop_mps
  ci=0; : > "$raw"
  for cfg in "${sel[@]}"; do
    ci=$((ci+1))
    awk -v c=$ci '/^PROF/{print c, $2, $3, $4, $5} /^PION/{print c, "PION", $2, $3}' "$cfg.bag.$tag" >> "$raw"
  done
  py_args+=("$m=$raw")
done

# --- analyze: per-mass bag s_T and m_pi, then the chiral extrapolation -----------------------------
echo "== analyze: the bag shape s_T(m_pi) and the chiral trend =="
python3 - "$T" "$PLAT_LO" "$PLAT_HI" "$R0A" "$L" "$MAXR0FRAC" "$MIN_MPI" "${py_args[@]}" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
T, lo, hi, r0a = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
L, maxr0frac, min_mpi = int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
items = [a.split("=", 1) for a in sys.argv[8:]]          # [(mass, file), ...]
r0_cut = maxr0frac * L

pts, table = [], []
single = None
for mass, fn in items:
    rows = [l.split() for l in open(fn) if l.split()]
    prof = np.array([[float(x) for x in r] for r in rows if r[1] != "PION" and len(r) == 5])
    pion = np.array([[float(r[0]), float(r[2]), float(r[3])] for r in rows if r[1] == "PION"])
    res = lat.bag_profile(prof, T=T, plateau=(lo, hi), r0_over_a=r0a)
    mpi = lat.correlator_mass(pion, T=T, tmin=lo, tmax=hi)
    single = (res, mpi, mass)
    # quality flag: a delocalised bag (fills the box) or a near-massless pion = exceptional config
    ok = (res["R0"] < r0_cut) and (mpi["mass"] > min_mpi)
    table.append((float(mass), mpi["mass"], res["R0"], res["s_T"], res["s_T_err"], ok))
    if ok:
        pts.append((mpi["mass"] ** 2, res["s_T"], res["s_T_err"]))

print(f"  {'m_q':>7} {'m_pi a':>9} {'R0/a':>7} {'s_T=R0/r0':>13}  window[0.43,0.70]")
for mq, mp, R0, sT, sTe, ok in sorted(table, key=lambda z: z[1]):
    where = "IN" if 0.43 <= sT <= 0.70 else ("below" if sT < 0.43 else "above")
    flag = "" if ok else "  <- EXCEPTIONAL (delocalised / massless): DROPPED"
    print(f"  {mq:>7.3f} {mp:>9.4f} {R0:>7.2f} {sT:>7.3f}+/-{sTe:<5.3f}  {where}{flag}")
ndrop = sum(1 for *_, ok in table if not ok)
if ndrop:
    print(f"  ({ndrop} exceptional point(s) excluded from the chiral fit: R0 > {r0_cut:.1f}a or "
          f"m_pi a < {min_mpi})")

if len(pts) >= 2:
    tr = lat.bag_chiral_trend(pts, r0_over_a=r0a)
    print(f"\n  chiral extrapolation s_T = c0 + c1 m_pi^2 :  slope c1 = {tr['slope']:+.3f}, "
          f"rises toward chiral = {tr['rising']}")
    print(f"  chiral s_T(m_pi=0) = {tr['chiral_s_T']:.3f} +/- {tr['chiral_s_T_err']:.3f} r0   "
          f"(lightest measured s_T = {tr['lightest_s_T']:.3f}, in-window={tr['any_in_window']})")
    print(f"  VERDICT: {tr['verdict']}")
elif len(pts) == 1:
    print(f"\n  only one clean mass survived the quality cut -- need >=2 to extrapolate; "
          f"loosen the scan or the cut.")
else:
    res, mpi, mass = single
    print(f"\n  single mass {mass} (m_pi a={mpi['mass']:.3f}): {res['verdict']}")
PY
