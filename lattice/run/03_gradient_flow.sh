#!/usr/bin/env bash
# Target 3: scale setting via the gradient (Wilson) flow -> w0 (cheap, minutes).
# Flow each config, record t^2<E>(t); w0^2 is where W(t)=t d/dt[t^2 E]=0.3.
set -euo pipefail
source "$(dirname "$0")/_lib.sh"                    # sets GRID, require_exe, note_params
OUT="${OUT:-out/flow}"; mkdir -p "$OUT"
ENS="${ENS:-out/puregauge}"   # reuse the target-1 ensemble (or any ensemble)
require_exe "$GRID/tests/smearing/Test_WilsonFlow"

L=24 ; T=48
: > "$OUT/t2E_raw.dat"        # columns: t  t^2<E>  (per config)
for cfg in "$ENS"/cfg.*; do
  "$GRID/tests/smearing/Test_WilsonFlow" --config "$cfg" --grid $L.$L.$L.$T \
      --flowtime 4.0 --flowsteps 400 >> "$OUT/t2E_raw.dat"
done
python3 - "$OUT/t2E_raw.dat" "$OUT/t2E.dat" <<'PY'
import sys, numpy as np
raw=np.loadtxt(sys.argv[1]); t=np.unique(raw[:,0])
np.savetxt(sys.argv[2], [[tt, raw[raw[:,0]==tt,1].mean()] for tt in t], header="t  t2E")
PY

echo "== analyze: w0 scale =="
python3 - "$OUT/t2E.dat" <<'PY'
import sys, numpy as np
from cartasis_sims import lattice as lat
t, t2E = np.loadtxt(sys.argv[1], unpack=True)
w0 = lat.gradient_flow_w0(t, t2E, ref=0.3)
print(f"  w0 = {w0:.4f} a   ->  a = w0_phys / w0   (w0_phys ~ 0.17 fm convention)")
PY
