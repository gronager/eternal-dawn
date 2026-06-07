#!/usr/bin/env bash
# Target L4 (production): the torsiton with DYNAMICAL fundamental sea quarks. The quenched pilot
# (run/06) proved the torsiton binds and showed the chiral trend, but on a frozen vacuum. This run
# includes the fermion determinant -- the real, fluctuating, chiral-symmetry-breaking condensate --
# so the torsiton mass, the m_N/sqrt(sigma) scale, and any excited rungs are the physical numbers.
#
# (1) generate a dynamical SU(3) N_f=2 fundamental ensemble (generate_dynamical), (2) measure the
# pion + nucleon (the torsiton) on it (the SAME run/06 spectroscopy, which is sea-agnostic), (3) read
# the masses off the plateaus. Cost: this is the leadership-class step (the CG sea force makes each
# trajectory far dearer than pure gauge) -- start small, watch acceptance and dH, scale up.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, MPS, run_pool
HERE="$(cd "$(dirname "$0")" && pwd)"
GEN="$HERE/../src/generate_dynamical"
require_exe "$GEN"

L="${L:-16}" ; T="${T:-32}" ; GRIDSPEC="$L.$L.$L.$T"
BETA="${BETA:-5.6}"
MASS="${MASS:--0.5}"        # dynamical sea-quark mass (lighter than quenched m_c: the sea suppresses
                           # exceptional configs). Tune toward the dynamical critical point.
NTRAJ="${NTRAJ:-1000}"      # total trajectories (dynamical thermalises slower than pure gauge)
MDSTEPS="${MDSTEPS:-20}"    # raise if acceptance is low (the fermion force is stiff)
SAVE="${SAVE:-5}"          # checkpoint cadence
NSTREAMS="${NSTREAMS:-1}"  # independent streams (raise to fill the GPU once one stream is tuned)
OUT="${OUT:-out/dyn_L${L}x${T}_m${MASS}}"; mkdir -p "$OUT"

# --- 1. generate the dynamical ensemble (skip if configs already exist) ---------------------
shopt -s nullglob
_existing=( "$OUT"/ckpoint_lat.* "$OUT"/stream*/ckpoint_lat.* )
if [ "${#_existing[@]}" -gt 0 ]; then
  echo "using existing dynamical configs in $OUT (${#_existing[@]} found)"
else
  echo "generating $NSTREAMS x $NTRAJ dynamical trajectories (beta=$BETA mass=$MASS $GRIDSPEC mdsteps=$MDSTEPS)"
  echo "  WATCH: acceptance (~0.7-0.9 is healthy; raise MDSTEPS if low) and dH per trajectory (~O(1))."
  start_mps
  for k in $(seq 1 "$NSTREAMS"); do
    d="$OUT/stream$k"
    printf 'mkdir -p %q && cd %q && %q --grid %q --beta %q --mass %q --seed %d --mdsteps %d --save-interval %d --StartingType HotStart --Trajectories %d --accelerator-threads 8 > hmc.log 2>&1\n' \
      "$d" "$d" "$GEN" "$GRIDSPEC" "$BETA" "$MASS" "$((k*1000 + 1))" "$MDSTEPS" "$SAVE" "$NTRAJ"
  done | run_pool "$NSTREAMS"
  stop_mps
  n=$(ls "$OUT"/stream*/ckpoint_lat.* 2>/dev/null | grep -vE '\.gz$' | wc -l)
  echo "done: $n configs in $OUT/stream*/"
fi

# --- 2. measure the torsiton on the dynamical ensemble (same spectroscopy as run/06) --------
# Run the valence at the SAME mass as the sea for the unitary point (or scan VALMASS for a partially-
# quenched chiral fit). The run/06 analysis (lattice.baryon_spectrum) reads the plateaus.
VALMASS="${VALMASS:-$MASS}"
echo "== measuring the torsiton (valence mass $VALMASS) on the dynamical ensemble =="
echo "   ->  OUT=$OUT L=$L T=$T MASS=$VALMASS THERM=<therm> STRIDE=<stride> ./run/06_baryon_spectrum.sh"
echo "   (set THERM past the dynamical thermalisation -- check stream*/hmc.log for the plaquette"
echo "    settling -- and STRIDE past the (longer) dynamical autocorrelation.)"
