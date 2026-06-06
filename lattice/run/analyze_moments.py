#!/usr/bin/env python3
"""Turn measure_modenumber's Chebyshev-moment output into nu(M) and gamma_m.

  analyze_moments.py moments.txt --mlo 0.20 --mhi 0.50 [--free L T]

Parses the 'XMAX <v>' and 'MOMENT <k> <mu>' lines, builds nu(M) via the validated
mode_number_from_chebyshev_moments (Jackson-damped KPM), and fits gamma_m over [mlo, mhi].
With --free it also prints the analytic free Wilson gamma_m and the nu(M) agreement on the
same window -- the known-answer check for a --free measurement run.
"""
import argparse
import sys

import numpy as np

from cartasis_sims import lattice as lat


def parse_moments(path):
    xmax, mu = None, {}
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) >= 2 and tok[0] == "XMAX":
                xmax = float(tok[1])
            elif len(tok) >= 3 and tok[0] == "MOMENT":
                mu[int(tok[1])] = float(tok[2])
    if xmax is None or not mu:
        sys.exit("no XMAX / MOMENT lines found -- is this measure_modenumber output?")
    return xmax, np.array([mu[k] for k in range(max(mu) + 1)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--mlo", type=float, default=0.20)
    ap.add_argument("--mhi", type=float, default=0.50)
    ap.add_argument("--npts", type=int, default=40)
    ap.add_argument("--free", nargs=2, type=int, metavar=("L", "T"),
                    help="also print the analytic free Wilson answer on an L^3 x T lattice")
    a = ap.parse_args()

    xmax, mu = parse_moments(a.file)
    M = np.linspace(a.mlo * 0.8, a.mhi * 1.2, a.npts)
    nu = lat.mode_number_from_chebyshev_moments(mu, M, xmax)
    out = lat.anomalous_dimension_from_mode_number(M, nu, window=(a.mlo, a.mhi))
    print(f"xmax={xmax:.4f}  moments={len(mu) - 1}  mu0={mu[0]:.0f} (= 12*sites)")
    print(f"gamma_m = {out['gamma_m']:+.4f}  (slope {out['slope']:.3f})  on M in [{a.mlo}, {a.mhi}]")

    if a.free:
        nuf = lat.free_wilson_mode_number(a.free[0], a.free[1], M)
        of = lat.anomalous_dimension_from_mode_number(M, nuf, window=(a.mlo, a.mhi))
        rel = np.abs(nu - nuf) / np.maximum(nuf, 1)
        print(f"free analytic gamma_m = {of['gamma_m']:+.4f}  (KPM should match this)")
        print(f"nu(M) KPM vs analytic free: median rel.err = {np.median(rel):.3f}")


if __name__ == "__main__":
    main()
