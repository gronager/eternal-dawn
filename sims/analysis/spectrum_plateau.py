#!/usr/bin/env python3
r"""Plateau spectroscopy for the L4 ensemble: correlators -> m_pi a, m_N a, m_N/m_pi.

Reads two-point correlator files (pion and/or nucleon), builds jackknife effective
masses, fits a constant plateau over a chosen window, and reports the masses and --
the number the bag sharpness s_T keys on -- the ratio m_N/m_pi with a correlated
jackknife error. It deliberately stops at the lattice observables; feed the printed
ratio into the existing run/10 ground-state line for s_T itself.

You do NOT need to wait for trajectory 1500: point this at correlators measured on the
thermalised configs you already have. Ratios converge faster than absolute masses, so a
few tens of configs already bracket m_N/m_pi.

Input formats (auto-detected; all whitespace/CSV text)
------------------------------------------------------
  (a) cfg t  C            three columns, one row per (config, timeslice)
  (b) t  C                two columns; repeated blocks (t resets to 0) = configs
  (c) matrix              rows = configs, cols = timeslices (Nt columns)
A 4th column (imaginary part) is ignored; the real part (col 'C') is used.

Usage
-----
  python spectrum_plateau.py --pion pion.dat --nt 64 --type cosh --tmin 12 --tmax 24
  python spectrum_plateau.py --pion pion.dat --nucleon nucl.dat --nt 64 \
        --pion-type cosh --nucleon-type exp --pion-window 12 24 --nucleon-window 8 16
  python spectrum_plateau.py --pion pion.dat            # auto window, prints effmass table

If you omit --tmin/--tmax (or --*-window) the script auto-picks the flattest plateau and
tells you what it chose -- then you pin it explicitly for the final number.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np


# ---- loading --------------------------------------------------------------------------
def load_corr(path, nt=None):
    """Return C[config, t] (real). Auto-detects the three layouts above."""
    a = np.loadtxt(path, ndmin=2)
    ncol = a.shape[1]
    if ncol >= 3 and _looks_like_index(a[:, 0]):
        # (cfg, t, C[, im])
        cfgs = np.unique(a[:, 0])
        ts = np.unique(a[:, 1])
        nt_ = len(ts)
        C = np.full((len(cfgs), nt_), np.nan)
        tindex = {t: i for i, t in enumerate(sorted(ts))}
        cindex = {c: i for i, c in enumerate(sorted(cfgs))}
        for row in a:
            C[cindex[row[0]], tindex[row[1]]] = row[2]
        return C
    if ncol == 2 or (ncol >= 2 and not _looks_like_index(a[:, 0])):
        # (t, C) repeated blocks: split where t resets
        t = a[:, 0]
        starts = np.where(t == t[0])[0]
        if len(starts) > 1:
            blocks = np.split(a, starts[1:])
            nt_ = len(blocks[0])
            C = np.array([b[:nt_, 1] for b in blocks if len(b) >= nt_])
            return C
        return a[:, 1][None, :]  # single config (no jackknife possible)
    # matrix: rows=configs, cols=t  (optionally trim to nt)
    return a[:, :nt] if nt else a


def _looks_like_index(col):
    """Heuristic: a config/time index column is non-negative integers with repeats."""
    v = np.asarray(col)
    if not np.allclose(v, np.round(v)):
        return False
    return (v.min() >= 0) and (len(np.unique(v)) < len(v))


# ---- jackknife ------------------------------------------------------------------------
def jack_samples(C):
    """Leave-one-out jackknife resamples of the correlator: returns Cjk[i, t]."""
    n = C.shape[0]
    tot = C.sum(0)
    return (tot[None, :] - C) / (n - 1)


def eff_mass(Cmean, kind="cosh", nt=None):
    """Effective mass array m_eff(t) for a mean correlator."""
    Cm = np.asarray(Cmean, float)
    T = len(Cm)
    with np.errstate(all="ignore"):
        if kind == "cosh":
            m = np.full(T, np.nan)
            for t in range(1, T - 1):
                r = (Cm[t - 1] + Cm[t + 1]) / (2 * Cm[t])
                m[t] = np.arccosh(r) if r >= 1 else np.nan
            return m
        # exponential (forward): log C(t)/C(t+1)
        m = np.full(T, np.nan)
        m[:-1] = np.log(Cm[:-1] / Cm[1:])
        return m


def eff_mass_jack(C, kind="cosh"):
    """Jackknife mean & error of the effective-mass curve."""
    Cjk = jack_samples(C)
    n = Cjk.shape[0]
    curves = np.array([eff_mass(Cjk[i], kind) for i in range(n)])
    mean = np.nanmean(curves, axis=0)
    err = np.sqrt((n - 1) / n * np.nansum((curves - mean) ** 2, axis=0))
    return mean, err, curves


def plateau_fit(curves, tmin, tmax):
    """Constant (uncorrelated, inverse-variance) fit over [tmin, tmax] on jackknife curves.

    Returns (mass, err, chi2_dof). Operates per jackknife sample so the error is a proper
    jackknife error; chi2/dof uses the across-sample variance per timeslice."""
    n = curves.shape[0]
    sl = slice(tmin, tmax + 1)
    seg = curves[:, sl]                      # [n_jack, n_t]
    mean_t = np.nanmean(seg, axis=0)
    var_t = (n - 1) / n * np.nansum((seg - mean_t) ** 2, axis=0)
    var_t = np.where(var_t > 0, var_t, np.nan)
    w = 1.0 / var_t
    good = np.isfinite(w) & np.isfinite(mean_t)
    if good.sum() < 1:
        return np.nan, np.nan, np.nan
    # per-sample weighted-mean plateau value -> jackknife error
    persample = np.nansum(np.where(good, seg * w, 0.0), axis=1) / np.nansum(w[good])
    mass = persample.mean()
    err = np.sqrt((n - 1) / n * np.sum((persample - mass) ** 2))
    chi2 = np.nansum(((mean_t - mass) ** 2 * w)[good])
    dof = max(1, good.sum() - 1)
    return mass, err, chi2 / dof


def auto_window(mean, err, kind="cosh", minlen=4):
    """Pick the flattest run of timeslices (smallest spread / largest length)."""
    T = len(mean)
    lo = 2 if kind == "cosh" else 1
    hi = T - 2
    valid = np.where(np.isfinite(mean[lo:hi]) & np.isfinite(err[lo:hi]))[0] + lo
    if len(valid) < minlen:
        return (lo, min(lo + minlen, hi))
    best, best_score = (valid[0], valid[-1]), np.inf
    for i in range(len(valid)):
        for j in range(i + minlen - 1, len(valid)):
            a, b = valid[i], valid[j]
            seg = mean[a:b + 1]
            e = err[a:b + 1]
            if not np.all(np.isfinite(seg)):
                continue
            spread = np.nanstd(seg) / (1 + b - a)        # reward flat & long
            score = spread + np.nanmean(e) / np.sqrt(b - a + 1)
            if score < best_score:
                best_score, best = score, (a, b)
    return best


def analyse(name, path, nt, kind, window):
    C = load_corr(path, nt)
    if C.shape[0] < 2:
        print(f"  [{name}] only {C.shape[0]} config(s) -- need >=2 for jackknife.")
        return None
    mean, err, curves = eff_mass_jack(C, kind)
    if window is None:
        tmin, tmax = auto_window(mean, err, kind)
        chosen = "auto"
    else:
        tmin, tmax = window
        chosen = "manual"
    mass, merr, chi2dof = plateau_fit(curves, tmin, tmax)
    print(f"\n[{name}]  ({C.shape[0]} configs, Nt={C.shape[1]}, {kind})")
    print(f"  effective mass m_eff(t):")
    for t in range(len(mean)):
        if np.isfinite(mean[t]):
            mark = "  <" if (tmin <= t <= tmax) else ""
            print(f"    t={t:3d}   {mean[t]:8.4f} +/- {err[t]:7.4f}{mark}")
    print(f"  plateau [{tmin},{tmax}] ({chosen}):  "
          f"m a = {mass:.4f} +/- {merr:.4f}   chi2/dof = {chi2dof:.2f}")
    return dict(name=name, C=C, mean=mean, err=err, curves=curves,
                tmin=tmin, tmax=tmax, mass=mass, merr=merr, kind=kind)


def ratio_jack(num, den):
    """Correlated jackknife of m_N/m_pi (same configs assumed, aligned by index)."""
    n = min(num["curves"].shape[0], den["curves"].shape[0])
    # per-sample plateau values
    def per(d):
        sl = slice(d["tmin"], d["tmax"] + 1)
        seg = d["curves"][:n, sl]
        m = np.nanmean(seg, axis=0)
        var = np.nanvar(seg, axis=0)
        w = 1.0 / np.where(var > 0, var, np.nan)
        good = np.isfinite(w)
        return np.nansum(np.where(good, seg * w, 0.0), axis=1) / np.nansum(w[good])
    rN, rpi = per(num), per(den)
    r = rN / rpi
    mean = r.mean()
    err = np.sqrt((n - 1) / n * np.sum((r - mean) ** 2))
    return mean, err


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pion", help="pion correlator file")
    p.add_argument("--nucleon", help="nucleon correlator file")
    p.add_argument("--nt", type=int, default=None, help="time extent (for matrix input / cosh)")
    p.add_argument("--type", default="cosh", choices=["cosh", "exp"],
                   help="default fit form for both channels")
    p.add_argument("--pion-type", choices=["cosh", "exp"])
    p.add_argument("--nucleon-type", choices=["cosh", "exp"])
    p.add_argument("--tmin", type=int, help="plateau tmin (applies to both unless --*-window)")
    p.add_argument("--tmax", type=int)
    p.add_argument("--pion-window", type=int, nargs=2, metavar=("TMIN", "TMAX"))
    p.add_argument("--nucleon-window", type=int, nargs=2, metavar=("TMIN", "TMAX"))
    p.add_argument("--out", default="sims/output/lattice")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    if not (args.pion or args.nucleon):
        p.error("give --pion and/or --nucleon")

    base_win = (args.tmin, args.tmax) if (args.tmin is not None and args.tmax is not None) else None
    results = {}
    if args.pion:
        win = tuple(args.pion_window) if args.pion_window else base_win
        results["pion"] = analyse("pion", args.pion, args.nt,
                                  args.pion_type or args.type, win)
    if args.nucleon:
        win = tuple(args.nucleon_window) if args.nucleon_window else base_win
        results["nucleon"] = analyse("nucleon", args.nucleon, args.nt,
                                     args.nucleon_type or args.type, win)

    print("\n=== summary ====================================================")
    pi = results.get("pion")
    nu = results.get("nucleon")
    if pi:
        print(f"  m_pi a = {pi['mass']:.4f} +/- {pi['merr']:.4f}")
    if nu:
        print(f"  m_N  a = {nu['mass']:.4f} +/- {nu['merr']:.4f}")
    if pi and nu:
        rm, re_ = ratio_jack(nu, pi)
        print(f"  m_N / m_pi = {rm:.4f} +/- {re_:.4f}   "
              f"(constituent-counting target 3/2 = 1.5)")
        print("  -> feed this ratio into the run/10 ground-state line for s_T.")
    print("================================================================")

    if not args.no_plots:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"(plots skipped: {e})")
            return 0
        os.makedirs(args.out, exist_ok=True)
        chans = [r for r in (pi, nu) if r]
        if chans:
            fig, axes = plt.subplots(1, len(chans), figsize=(5.2 * len(chans), 4),
                                     squeeze=False)
            for ax, d in zip(axes[0], chans):
                t = np.arange(len(d["mean"]))
                ax.errorbar(t, d["mean"], yerr=d["err"], fmt="o", ms=3, color="C0")
                ax.axhspan(d["mass"] - d["merr"], d["mass"] + d["merr"],
                           xmin=0, xmax=1, color="C1", alpha=0.3)
                ax.axhline(d["mass"], color="C1",
                           label=f"m a = {d['mass']:.3f}({d['merr']:.3f})")
                ax.axvspan(d["tmin"], d["tmax"], color="0.85", alpha=0.5)
                ax.set_title(f"{d['name']} effective mass ({d['kind']})")
                ax.set_xlabel("t"); ax.set_ylabel(r"$m_{\rm eff}\,a$")
                ax.legend(fontsize=8); ax.grid(alpha=0.2)
                finite = d["mean"][np.isfinite(d["mean"])]
                if finite.size:
                    ax.set_ylim(max(0, finite.min() * 0.5), finite.max() * 1.5)
            fig.tight_layout()
            out = os.path.join(args.out, "spectrum_plateau.png")
            fig.savefig(out, dpi=130); plt.close(fig)
            print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
