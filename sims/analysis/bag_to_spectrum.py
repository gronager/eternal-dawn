#!/usr/bin/env python3
r"""Measured torsiton bag -> generation spectrum (the honest, no-sharp-source-assumption read).

Reads a bag-profile file (columns `r/a  rho(r)`, the output of run/09_bag_profile or
run/10_condensate_3pt, e.g. figures/data/torsiton_bag_L16x48.dat), and runs the corrected
analysis built in cartasis_sims.well_spectrum:

  1. fit the SHAPE in full -- not just the half-density radius R0, but the WALL thickness a_wall
     and the tail length, because (well_spectrum, fixed-(R0,depth) test) the generation SPAN is
     carried by the wall, not the radius: a soft wall gives a small span;
  2. build the confining well  M(r) = m_vac * (1 - rho(r))  (Picture A: chiral-restored core,
     vacuum mass outside) on a relativistic grid -- the FULL measured shape, no Woods-Saxon fit;
  3. solve the first-order (G,F) Dirac (doubler-free matched shooting) for the bound spectrum;
  4. report the level COUNT (set by the depth m_vac*r0) and the configurational-mass SPAN (set by
     the wall), against the observed lepton span ~3477.

The point: this solves the Dirac equation in the MEASURED bag. It does NOT assume a sharp Gaussian
source (the overlap model that the chapter's headline span relies on). Where the two disagree is
exactly the open dictionary between the measured condensate profile and the effective scalar mass.

Usage
-----
  python bag_to_spectrum.py ../figures/data/torsiton_bag_L16x48.dat --r0a 3.131 --mvac-r0 1.47
  # depth from the nucleon mass instead (m_vac*r0 = (m_N a / 3) * r0/a):
  python bag_to_spectrum.py bag.dat --r0a 3.131 --mN 1.4
  # scan the depth to see what 3 levels would require:
  python bag_to_spectrum.py bag.dat --r0a 3.131 --scan 1.5 2 4 8
"""
from __future__ import annotations

import argparse
import numpy as np

from cartasis_sims import well_spectrum as ws

LEPTON_SPAN = 3477.0


def load_bag(path):
    """Load (r/a, rho) and any `r0/a = ...` hint from the header comments."""
    r, rho, r0a = [], [], None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if "r0/a" in s:
                after = s.split("r0/a", 1)[1].lstrip(" =\t")
                for tok in after.replace(",", " ").split():
                    tok = tok.strip(".;)")             # trailing punctuation, e.g. "3.131."
                    try:
                        r0a = float(tok); break
                    except ValueError:
                        continue
            continue
        a, b = s.split()[:2]
        r.append(float(a)); rho.append(float(b))
    r, rho = np.asarray(r), np.asarray(rho)
    if rho[0] != 0:
        rho = rho / rho[0]
    return r, rho, r0a


def fit_shape(r, rho):
    """Full shape: Fermi wall (R0, a_wall) and the large-r exponential tail length lam (lattice a)."""
    from scipy.optimize import curve_fit
    fermi = lambda rr, R0, a: 1.0 / (1.0 + np.exp((rr - R0) / a))
    (R0, a_wall), _ = curve_fit(fermi, r, rho, p0=[1.8, 0.6], maxfev=40000)
    big = r > max(2.0, 0.5 * r.max())
    lam = float("nan")
    if big.sum() >= 3 and np.all(rho[big] > 0):
        slope = np.polyfit(r[big], np.log(rho[big]), 1)[0]
        lam = -1.0 / slope if slope < 0 else float("nan")
    return float(R0), float(a_wall), lam


def spectrum_at_depth(r_r0, rho, mvac_r0, N=900, R=12.0):
    """Build M(r)=mvac_r0*(1-rho) on a solver grid (lengths in r0) and solve the Dirac spectrum."""
    rg = np.linspace(R / N, R, N)
    M = np.interp(rg, r_r0, mvac_r0 * (1.0 - rho), left=0.0, right=mvac_r0)
    res = ws.spectrum_from_well(rg, M, np.zeros_like(rg))
    sl, sc = ws.spans(res)
    return res, sl, sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bagfile")
    ap.add_argument("--r0a", type=float, default=None, help="r0/a (else read from header)")
    ap.add_argument("--mvac-r0", type=float, default=None, help="depth m_vac*r0 directly")
    ap.add_argument("--mN", type=float, default=None, help="nucleon m_N a -> depth (m_N/3)*r0a")
    ap.add_argument("--scan", type=float, nargs="*", default=None, help="scan these m_vac*r0 depths")
    args = ap.parse_args()

    r_a, rho, r0a_hdr = load_bag(args.bagfile)
    r0a = args.r0a or r0a_hdr
    if r0a is None:
        ap.error("need r0/a (--r0a or a header line)")
    r_r0 = r_a / r0a

    R0, a_wall, lam = fit_shape(r_a, rho)
    s_T = R0 / r0a
    print(f"bag: {len(r_a)} points, r0/a={r0a}")
    print(f"  shape:  R0={R0:.2f}a  s_T=R0/r0={s_T:.3f}   wall a_wall={a_wall:.2f}a "
          f"(a_wall/R0={a_wall/R0:.2f})   tail lam={lam:.2f}a")
    soft = a_wall / R0 > 0.30
    print(f"  -> wall is {'SOFT (small span expected)' if soft else 'SHARP (large span possible)'} "
          f"(the span lever is the wall, not R0)")

    depths = []
    if args.mvac_r0 is not None:
        depths.append((args.mvac_r0, "given"))
    if args.mN is not None:
        depths.append(((args.mN / 3.0) * r0a, f"m_N/3 (m_N a={args.mN})"))
    for d in (args.scan or []):
        depths.append((d, "scan"))
    if not depths:
        depths = [(1.47, "default")]

    print(f"\n  {'m_vac*r0':>9}  {'n_bound':>7}  {'span_local':>10}  {'span_core':>9}   E_n/M_vac   [source]")
    for mvac_r0, src in depths:
        res, sl, sc = spectrum_at_depth(r_r0, rho, mvac_r0)
        es = " ".join(f"{x['E']/res['M_vac']:.3f}" for x in res["levels"])
        sls = f"{sl:10.1f}" if sl else f"{'-':>10}"
        scs = f"{sc:9.1f}" if sc else f"{'-':>9}"
        print(f"  {mvac_r0:>9.2f}  {res['n_bound']:>7}  {sls}  {scs}   {es}   [{src}]")
    print(f"\n  observed lepton span = {LEPTON_SPAN:.0f}; the consecutive ladder needs ~3 bound levels.")
    print("  NB: this is the Dirac spectrum IN the measured well (configurational mass), with NO sharp-")
    print("  source assumption. A larger span requires a sharper measured wall or the condensate->mass")
    print("  dictionary to be sharper than M=m_vac(1-rho).")


if __name__ == "__main__":
    main()
