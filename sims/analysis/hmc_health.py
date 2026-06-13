#!/usr/bin/env python3
r"""HMC health check for the L4 light-sea ensemble (Grid HMC log -> verdict + plots).

Tuned to Grid's actual HMC summary format (GridLogHMC channel):

    Total H before trajectory = <H0>
    Total H after trajectory  = <H1>  dH = <dH>
    Skipping Metropolis test                         <- thermalisation (forced accept)
    exp(-dH) = <prob>   Random = <rn>                <- production
    Metropolis_test -- ACCEPTED | REJECTED

Plaquette is NOT logged by name, but the Wilson gauge action is
(`S [1][0] H = <Sg>`), so we reconstruct it:  <plaq> = 1 - Sg/(beta * 6 * V).
beta and the lattice volume are auto-read from the log banner
("beta=5.6 ..." and "Full Dimensions : [16 16 16 64]"); override with --beta/--volume.

Reports
-------
  * thermalisation: how many trajectories Grid skipped Metropolis on, and where the
    plaquette plateau actually sets in (the cut that matters for measurements)
  * acceptance over the PRODUCTION trajectories (want ~0.70-0.90)
  * <exp(-dH)> == 1 over production (the unbiasedness check)
  * <dH>, dH distribution
  * plaquette plateau mean, first/second-half drift, integrated autocorr time

Usage
-----
  python hmc_health.py lattice/out/dyn_L16x64_m-0.75/stream1/hmc.log
  python hmc_health.py hmc.log --therm 100          # pin the measurement cut (trajs)
  python hmc_health.py hmc.log --beta 5.6 --volume 262144   # if banner not present

Exit 0 if healthy, 1 if a WARN fires, 2 if the parse is too thin to judge.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

RE_BETA = re.compile(r"beta=([0-9.]+)")
RE_DIMS = re.compile(r"Full Dimensions\s*:\s*\[([0-9 ]+)\]")
RE_HBEFORE = re.compile(r"before trajectory\s*=\s*([-+0-9.eEdD]+)")
RE_HAFTER = re.compile(r"after trajectory\s*=\s*([-+0-9.eEdD]+)\s*dH\s*=\s*([-+0-9.eEdD]+)")
RE_EXPMDH = re.compile(r"exp\(-dH\)\s*=\s*([-+0-9.eEdD]+)")
RE_METRO = re.compile(r"Metropolis_test\s*--\s*(ACCEPTED|REJECTED)", re.IGNORECASE)
RE_GAUGE = re.compile(r"S\s*\[1\]\[0\]\s*H\s*=\s*([-+0-9.eEdD]+)")


def _f(s):
    return float(s.replace("d", "e").replace("D", "e"))


def parse_grid_log(path):
    recs, cur, gauge = [], {}, []
    beta, dims = None, None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if beta is None:
                m = RE_BETA.search(line)
                if m:
                    beta = float(m.group(1))
            if dims is None:
                m = RE_DIMS.search(line)
                if m:
                    dims = [int(x) for x in m.group(1).split()]

            if "before trajectory" in line:
                m = RE_HBEFORE.search(line)
                if m:
                    if cur:
                        recs.append(cur)
                    cur = {"H0": _f(m.group(1))}
            elif "after trajectory" in line:
                m = RE_HAFTER.search(line)
                if m:
                    cur["H1"] = _f(m.group(1))
                    cur["dH"] = _f(m.group(2))
            elif "Skipping Metropolis" in line:
                cur["acc"] = 1.0
                cur["therm"] = True
                recs.append(cur)
                cur = {}
            elif "Metropolis_test" in line:
                m = RE_METRO.search(line)
                if m:
                    cur["acc"] = 1.0 if m.group(1).upper() == "ACCEPTED" else 0.0
                cur["therm"] = False
                recs.append(cur)
                cur = {}
            else:
                m = RE_EXPMDH.search(line)
                if m:
                    cur["expmdH"] = _f(m.group(1))
            g = RE_GAUGE.search(line)
            if g:
                gauge.append(_f(g.group(1)))
    if cur:
        recs.append(cur)
    vol = int(np.prod(dims)) if dims else None
    return recs, np.array(gauge, float), beta, vol, dims


def jackknife(x, func=np.mean):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return (func(x) if n else np.nan), np.nan
    drop = np.array([func(np.delete(x, i)) for i in range(n)])
    return func(x), np.sqrt((n - 1) / n * np.sum((drop - drop.mean()) ** 2))


def tau_int(series, c=5.0):
    x = np.asarray(series, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 8:
        return np.nan, np.nan
    x = x - x.mean()
    var = np.dot(x, x) / n
    if var == 0:
        return 0.5, 0.0
    rho = np.array([np.dot(x[: n - t], x[t:]) / (n * var) for t in range(n)])
    tau, W = 0.5, 1
    for W in range(1, n):
        tau = 0.5 + np.sum(rho[1:W + 1])
        if W >= c * tau:
            break
    return tau, tau * np.sqrt((2.0 * (2 * W + 1)) / n)


def auto_plateau_cut(plaq):
    """First index from which the running mean STAYS within 1 sd of the 2nd-half mean."""
    if plaq.size < 20:
        return 0
    run = np.cumsum(plaq) / np.arange(1, plaq.size + 1)
    tail = plaq[plaq.size // 2:]
    band = max(tail.std(ddof=1), 1e-9)
    inside = np.abs(run - tail.mean()) <= band
    for i in range(plaq.size):
        if inside[i:].all():
            return i
    return plaq.size // 5


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("log")
    p.add_argument("--beta", type=float, default=None, help="override gauge beta")
    p.add_argument("--volume", type=int, default=None, help="override lattice volume (sites)")
    p.add_argument("--therm", type=int, default=None,
                   help="measurement thermalisation cut in TRAJECTORIES (default: auto plateau)")
    p.add_argument("--out", default="sims/output/lattice")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    recs, gauge, beta, vol, dims = parse_grid_log(args.log)
    beta = args.beta or beta
    vol = args.volume or vol
    n = len(recs)

    dH = np.array([r.get("dH", np.nan) for r in recs], float)
    acc = np.array([r.get("acc", np.nan) for r in recs], float)
    therm = np.array([r.get("therm", False) for r in recs], bool)
    n_dh = np.isfinite(dH).sum()

    print(f"\nParsed {n} trajectory records from {args.log}")
    print(f"  dH={n_dh}  acc={np.isfinite(acc).sum()}  gauge-action samples={gauge.size}  "
          f"beta={beta}  volume={vol} {dims if dims else ''}")

    if n < 5 or n_dh < 5:
        print("\nINSUFFICIENT DATA: parsed too few complete trajectories to judge. "
              "If the run has many trajectories, paste ~40 summary lines and I'll adjust.",
              file=sys.stderr)
        return 2

    # ---- plaquette from the gauge action -------------------------------------------------
    plaq = None
    if gauge.size and beta and vol:
        plaq_all = 1.0 - gauge / (beta * 6.0 * vol)
        # initial+final action per traj -> ~2 samples/traj; take one per traj
        plaq = plaq_all[::2] if gauge.size > 1.5 * n else plaq_all
    n_skip = int(therm.sum())

    print("\n--- HMC health -------------------------------------------------")
    print(f"[INFO] thermalisation: Grid skipped Metropolis on {n_skip} trajectories "
          f"(forced-accept hot-start relaxation)")

    warn = False
    prod = ~therm  # production = Metropolis applied

    # ---- acceptance (production only) ----------------------------------------------------
    ap = acc[prod & np.isfinite(acc)]
    if ap.size:
        rate = ap.mean()
        flag = "OK  " if 0.65 <= rate <= 0.92 else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[{flag}] acceptance (prod)  = {rate:6.3f}  over {ap.size} trajs   (want 0.70-0.90)")
    else:
        print("[INFO] acceptance        = (still in thermalisation; no Metropolis trajs yet)")

    # ---- <exp(-dH)> == 1 (production) ----------------------------------------------------
    dp = dH[prod & np.isfinite(dH)]
    if dp.size >= 5:
        em, ee = jackknife(np.exp(-dp))
        nsig = abs(em - 1.0) / ee if ee else np.nan
        flag = "OK  " if (np.isfinite(nsig) and nsig < 3) else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[{flag}] <exp(-dH)> (prod)  = {em:6.4f} +/- {ee:6.4f}  "
              f"({nsig:.1f}sigma from 1; want <3)")
        dm, de = jackknife(dp)
        print(f"       <dH> (prod)       = {dm:+6.4f} +/- {de:6.4f}   "
              f"(max |dH| = {np.abs(dp).max():.2f})")
    else:
        print("[INFO] <exp(-dH)>        = (need more production trajectories)")

    # ---- plaquette plateau / drift / autocorrelation -------------------------------------
    if plaq is not None and plaq.size > 10:
        cut = args.therm if args.therm is not None else auto_plateau_cut(plaq)
        cut = max(0, min(cut, plaq.size - 5))
        pl = plaq[cut:]
        pm, pe = jackknife(pl)
        h = pl.size // 2
        m1 = pl[:h].mean(); m2 = pl[h:].mean()
        e12 = np.hypot(pl[:h].std(ddof=1) / np.sqrt(h), pl[h:].std(ddof=1) / np.sqrt(pl.size - h))
        nsig = abs(m1 - m2) / e12 if e12 > 0 else np.nan
        tau, terr = tau_int(pl)
        flag = "OK  " if (np.isfinite(nsig) and nsig < 3) else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[INFO] plaquette plateau sets in ~traj {cut} "
              f"({'manual' if args.therm is not None else 'auto'})")
        print(f"[{flag}] plaquette (plateau)= {pm:8.6f} +/- {pe:8.6f}")
        print(f"       drift 1st vs 2nd half: {m1:.6f} vs {m2:.6f}  ({nsig:.1f}sigma; want <3)")
        if np.isfinite(tau):
            print(f"       tau_int(plaq)     = {tau:5.2f} +/- {terr:4.2f} traj  "
                  f"-> save every ~{max(1, int(round(2 * tau)))} traj for independence; "
                  f"N_eff ~ {pl.size / (2 * tau):.0f}")
    else:
        print("[INFO] plaquette         = (no gauge action / beta / volume to reconstruct)")

    print("----------------------------------------------------------------")
    print("VERDICT: WARN -- check the flagged line(s)." if warn
          else "VERDICT: all checks pass. Run looks healthy.")

    if not args.no_plots:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"(plots skipped: {e})")
            return 1 if warn else 0
        os.makedirs(args.out, exist_ok=True)
        panels = (["plaq"] if plaq is not None else []) + (["dH"] if n_dh > 0 else [])
        fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.7 * len(panels)),
                                 squeeze=False)
        axes = axes[:, 0]
        k = 0
        if plaq is not None:
            ax = axes[k]; k += 1
            ax.plot(np.arange(plaq.size), plaq, lw=0.7, color="C0")
            if plaq.size > 10:
                ax.axvline(cut, color="0.5", ls="--", label=f"plateau cut (~{cut})")
                ax.set_ylim(min(plaq.min(), pm - 6 * pe), pm + 6 * pe if pe else None)
                ax.legend(fontsize=8)
            ax.set_ylabel("plaquette\n(from gauge action)"); ax.grid(alpha=0.2)
            ax.set_title("HMC health: plaquette history (reconstructed)")
        if n_dh > 0:
            ax = axes[k]; k += 1
            ax.plot(np.arange(n), dH, lw=0.6, color="C3")
            ax.axhline(0, color="0.6", lw=0.6)
            ax.set_xlabel("trajectory"); ax.set_ylabel("dH"); ax.grid(alpha=0.2)
        fig.tight_layout()
        out = os.path.join(args.out, "hmc_health.png")
        fig.savefig(out, dpi=130); plt.close(fig)
        print(f"\nwrote {out}")

    return 1 if warn else 0


if __name__ == "__main__":
    raise SystemExit(main())
