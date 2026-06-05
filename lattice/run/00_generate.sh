#!/usr/bin/env bash
# Generate gauge configs as K INDEPENDENT HMC streams in parallel on one GPU (Eternal Dawn).
# Config generation is a Markov chain (can't parallelise within a stream), but independent
# streams -- different seeds, hot starts, own directories -- run side by side under MPS and
# pool their configs. This fills the GH200 during generation, where one serial chain leaves it
# idle. Uses our generate_gauge program (CLI beta + seed); build it with: make -C lattice/src
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # GRID, require_exe, start_mps/stop_mps, run_pool
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${OUT:-out/puregauge}"; mkdir -p "$OUT"
GEN="$HERE/../src/generate_gauge"
require_exe "$GEN"

L="${L:-16}" ; T="${T:-32}" ; GRIDSPEC="$L.$L.$L.$T"
BETA="${BETA:-5.6}"
NTRAJ="${NTRAJ:-200}"        # trajectories per stream
NSTREAMS="${NSTREAMS:-4}"    # independent parallel streams (raise to fill the GPU)

echo "generating $NSTREAMS independent streams x $NTRAJ trajectories (beta=$BETA, $GRIDSPEC)"
start_mps
for k in $(seq 1 "$NSTREAMS"); do
  d="$OUT/stream$k"
  # each stream: own dir, hot start, well-separated seed -> genuinely independent chain
  printf 'mkdir -p %q && cd %q && %q --grid %q --beta %q --seed %d --StartingType HotStart --Trajectories %d --accelerator-threads 8 > hmc.log 2>&1\n' \
    "$d" "$d" "$GEN" "$GRIDSPEC" "$BETA" "$((k*1000 + 1))" "$NTRAJ"
done | run_pool "$NSTREAMS"
stop_mps

n=$(ls "$OUT"/stream*/ckpoint_lat.* 2>/dev/null | grep -vE '\.gz$' | wc -l)
echo "done: $n configs across $OUT/stream*/  -- now run ./lattice/run/01_puregauge_potential.sh"
echo "(set NSTREAMS / NTRAJ / BETA via env; e.g. NSTREAMS=8 NTRAJ=400 ./lattice/run/00_generate.sh)"
