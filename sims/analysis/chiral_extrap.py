#!/usr/bin/env python3
r"""Chiral extrapolation of the torsiton (nucleon) from a valence-mass scan.

Reads several `baryon_raw_m-*.dat` files (columns: cfg t C_pi C_N -- the output of
run/06_baryon_spectrum.sh), fits m_pi a and m_N a per ensemble by jackknife plateau,
and extrapolates the trajectory:

  * GMOR line:      m_pi^2  vs  valence mass     ->  critical mass m_crit (m_pi^2 -> 0)
  * chiral nucleon: m_N     vs  m_pi^2           ->  m_N^0 (chiral-limit torsiton mass)

The pion is fit with the cosh effective mass (periodic, exact); the nucleon with the
forward log effective mass AND over several windows, because the window choice -- not the
statistics -- dominates the baryon error. Everything is in lattice units; multiply by 1/a
(from w0/r0, run/03) for MeV.

Usage
-----
  # masses inferred from the filenames (…m-0.75.dat -> -0.75):
  python chiral_extrap.py baryon_raw_m-0.75.dat baryon_raw_m-0.5.dat baryon_raw_m-0.3.dat

  # explicit masses / windows:
  python chiral_extrap.py f1 f2 f3 --masses -0.75 -0.5 -0.3 --T 64 \
        --pion-window 8 21 --nucleon-windows 10,20 12,24 16,26
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="invalid value encountered")


# ---- load + jackknife -----------------------------------------------------------------
def load_raw(path):
    """baryon_raw.dat (cfg t C_pi C_N) -> C_pi[cfg,t], C_N[cfg,t].

    Robust to non-numeric junk lines: Grid can leak a GPU-init line that starts with a
    digit (e.g. '0SharedMemoryNone:') past run/06's grep, so keep only rows of exactly four
    float-parseable columns."""
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) != 4:
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"{path}: no valid 'cfg t C_pi C_N' rows found")
    a = np.array(rows)
    cfgs = sorted(np.unique(a[:, 0]))
    nt = int(a[:, 1].max()) + 1
    ci = {c: i for i, c in enumerate(cfgs)}
    Cpi = np.full((len(cfgs), nt), np.nan)
    CN = np.full((len(cfgs), nt), np.nan)
    for r in a:
        i, t = ci[r[0]], int(r[1])
        Cpi[i, t], CN[i, t] = r[2], r[3]
    return Cpi, CN


def jack_samples(C):
    n = C.shape[0]
    return (C.sum(0)[None, :] - C) / (n - 1)


def eff_mass(Cm, kind):
    Cm = np.asarray(Cm, float)
    T = len(Cm)
    m = np.full(T, np.nan)
    with np.errstate(all="ignore"):
        if kind == "cosh":
            for t in range(1, T - 1):
                r = (Cm[t - 1] + Cm[t + 1]) / (2 * Cm[t])
                m[t] = np.arccosh(r) if r >= 1 else np.nan
        else:  # forward log (baryon): consistent-sign ratio is fine for negative C
            r = Cm[:-1] / Cm[1:]
            m[:-1] = np.where(r > 0, np.log(r), np.nan)
    return m


def eff_mass_jack(C, kind):
    Cjk = jack_samples(C)
    n = Cjk.shape[0]
    curves = np.array([eff_mass(Cjk[i], kind) for i in range(n)])
    mean = np.nanmean(curves, axis=0)
    err = np.sqrt((n - 1) / n * np.nansum((curves - mean) ** 2, axis=0))
    return mean, err, curves


def plateau(curves, tmin, tmax):
    """Inverse-variance constant fit over [tmin,tmax]; jackknife error + chi2/dof."""
    n = curves.shape[0]
    seg = curves[:, tmin:tmax + 1]
    mean_t = np.nanmean(seg, axis=0)
    var_t = (n - 1) / n * np.nansum((seg - mean_t) ** 2, axis=0)
    w = np.where(var_t > 0, 1.0 / var_t, np.nan)
    good = np.isfinite(w) & np.isfinite(mean_t)
    if good.sum() < 1:
        return np.nan, np.nan, np.nan
    persample = np.nansum(np.where(good, seg * w, 0.0), axis=1) / np.nansum(w[good])
    m = persample.mean()
    err = np.sqrt((n - 1) / n * np.sum((persample - m) ** 2))
    chi2 = np.nansum(((mean_t - m) ** 2 * w)[good]) / max(1, good.sum() - 1)
    return m, err, chi2


# ---- weighted straight-line fit (3 points, 1 dof) -------------------------------------
def wlinfit(x, y, dy):
    """y = a + b x, inverse-variance weighted. Returns (a, b, cov[2x2])."""
    x, y, dy = map(np.asarray, (x, y, dy))
    w = 1.0 / dy ** 2
    X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    beta = cov @ (X.T @ (w * y))
    return beta[0], beta[1], cov


def infer_mass(path):
    m = re.search(r"m(-?\d+(?:\.\d+)?)", os.path.basename(path))
    return float(m.group(1)) if m else None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="baryon_raw_m-*.dat files")
    p.add_argument("--masses", type=float, nargs="+",
                   help="valence masses (default: inferred from filenames)")
    p.add_argument("--T", type=int, default=None, help="time extent (default: from data)")
    p.add_argument("--pion-window", type=int, nargs=2, default=[8, 21], metavar=("TMIN", "TMAX"))
    p.add_argument("--nucleon-windows", nargs="+", default=["10,20", "12,24", "16,26"],
                   help="windows TMIN,TMAX for the nucleon systematic")
    p.add_argument("--out", default="sims/output/lattice")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    masses = args.masses or [infer_mass(f) for f in args.files]
    if any(m is None for m in masses):
        sys.exit("could not infer all masses from filenames; pass --masses")
    nwins = [tuple(int(x) for x in w.split(",")) for w in args.nucleon_windows]

    pts = []  # (mass, mpi, mpi_err, mN, mN_stat, mN_sys, ratio, ratio_err)
    print(f"{'mass':>7} {'m_pi a':>14} {'m_N a (stat) [sys]':>26} {'m_N/m_pi':>14}")
    print("-" * 70)
    for f, mv in zip(args.files, masses):
        Cpi, CN = load_raw(f)
        if Cpi.shape[0] < 2:
            print(f"{mv:>7.3f}  only {Cpi.shape[0]} config(s) -- skipped"); continue
        # pion (cosh)
        _, _, pic = eff_mass_jack(Cpi, "cosh")
        mpi, mpi_e, _ = plateau(pic, *args.pion_window)
        # nucleon (log) over several windows -> central + window systematic
        _, _, nuc = eff_mass_jack(CN, "log")
        mNs = [plateau(nuc, a, b) for a, b in nwins]
        vals = np.array([m for m, e, c in mNs if np.isfinite(m)])
        errs = np.array([e for m, e, c in mNs if np.isfinite(m)])
        mN = float(np.median(vals)); mN_stat = float(errs[np.argsort(vals)[len(vals)//2]])
        mN_sys = float((vals.max() - vals.min()) / 2) if len(vals) > 1 else 0.0
        ratio = mN / mpi
        ratio_e = ratio * np.hypot(mpi_e / mpi, np.hypot(mN_stat, mN_sys) / mN)
        pts.append((mv, mpi, mpi_e, mN, mN_stat, mN_sys, ratio, ratio_e))
        print(f"{mv:>7.3f}  {mpi:>7.4f}+/-{mpi_e:<6.4f}  "
              f"{mN:>7.4f}+/-{mN_stat:<6.4f}[+/-{mN_sys:.3f}]  {ratio:>6.3f}+/-{ratio_e:.3f}")

    if len(pts) < 2:
        sys.exit("\nneed >=2 mass points for the extrapolation")

    P = np.array(pts)
    mv, mpi, mpi_e, mN, mN_stat, mN_sys = (P[:, i] for i in range(6))
    mpi2 = mpi ** 2
    mpi2_e = 2 * mpi * mpi_e
    mN_tot = np.hypot(mN_stat, mN_sys)

    print("\n=== chiral extrapolation (lattice units) ======================")
    # GMOR: m_pi^2 = A (m - m_crit) -> fit m_pi^2 = a + b*m ; m_crit = -a/b
    a, b, cov = wlinfit(mv, mpi2, mpi2_e)
    m_crit = -a / b
    dmc = abs(m_crit) * np.hypot(np.sqrt(cov[0, 0]) / a, np.sqrt(cov[1, 1]) / b)
    print(f"GMOR  m_pi^2 = {b:.4f}*(m - m_crit),  m_crit = {m_crit:.4f} +/- {dmc:.4f}")
    print(f"      (your unitary sea mass {mv.min():.3f} sits "
          f"{(mv.min()-m_crit):.3f} above m_crit)")

    # chiral nucleon: m_N = m_N^0 + B*m_pi^2 -> intercept at m_pi^2=0
    a2, b2, cov2 = wlinfit(mpi2, mN, mN_tot)
    mN0, dmN0 = a2, np.sqrt(cov2[0, 0])
    print(f"m_N(m_pi^2=0) = m_N^0 = {mN0:.4f} +/- {dmN0:.4f}   "
          f"(slope dm_N/dm_pi^2 = {b2:.3f})")
    print("  -> multiply by 1/a (w0/r0, run/03) for MeV; this is the chiral torsiton mass.")
    print("================================================================")

    if not args.no_plots:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"(plots skipped: {e})"); return 0
        os.makedirs(args.out, exist_ok=True)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        xs = np.linspace(min(mv.min(), m_crit), mv.max(), 50)
        ax[0].errorbar(mv, mpi2, yerr=mpi2_e, fmt="o", color="C0")
        ax[0].plot(xs, a + b * xs, "C0-", lw=1)
        ax[0].axhline(0, color="0.6", lw=0.6); ax[0].axvline(m_crit, color="C3", ls="--",
                      label=f"$m_{{crit}}={m_crit:.3f}$")
        ax[0].set_xlabel("valence mass"); ax[0].set_ylabel(r"$(m_\pi a)^2$")
        ax[0].set_title("GMOR line"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2)
        xs2 = np.linspace(0, mpi2.max(), 50)
        ax[1].errorbar(mpi2, mN, yerr=mN_tot, fmt="o", color="C1")
        ax[1].plot(xs2, a2 + b2 * xs2, "C1-", lw=1)
        ax[1].errorbar([0], [mN0], yerr=[dmN0], fmt="*", ms=13, color="C3",
                       label=f"$m_N^0={mN0:.3f}({dmN0:.3f})$")
        ax[1].set_xlabel(r"$(m_\pi a)^2$"); ax[1].set_ylabel(r"$m_N a$")
        ax[1].set_title("chiral nucleon (torsiton)"); ax[1].legend(fontsize=8)
        ax[1].grid(alpha=0.2)
        fig.tight_layout()
        out = os.path.join(args.out, "chiral_extrap.png")
        fig.savefig(out, dpi=130); plt.close(fig)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
