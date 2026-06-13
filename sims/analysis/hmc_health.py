#!/usr/bin/env python3
r"""HMC health check for the L4 light-sea ensemble (Grid HMC log -> verdict + plots).

Reads a Grid HMC log (or a plain table) and reports the things that tell you, in one
glance, whether the run is sound:

  * acceptance rate                      (want ~0.70-0.90)
  * <exp(-dH)>  ==  1                     (THE unbiasedness check; must hold in expectation)
  * <dH>, dH distribution                (want O(1), no recurring spikes)
  * plaquette plateau + drift test       (thermalised?)
  * plaquette integrated autocorr time   (how many traj per independent config)

It is deliberately format-tolerant: it scans the log for the usual Grid lines
("Trajectory N", "Delta H = ...", "Average plaquette: ...", Metropolis Accept/Reject).
If your log differs, it still prints what it managed to parse so we can adjust the
regexes in one pass -- or use --table for a whitespace file with columns
(traj dH acc plaq).

Usage
-----
  python hmc_health.py RUN.log
  python hmc_health.py RUN.log --therm 200 --out sims/output/lattice
  python hmc_health.py ens.txt --table          # columns: traj dH acc plaq
  python hmc_health.py RUN.log --no-plots        # numbers only

Exit status is 0 if all checks pass, 1 if any WARN fires (handy in a cron).
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

# ---- default Grid log patterns (override-able if your build prints differently) ---------
RE_TRAJ = re.compile(r"Trajectory\s+(\d+)")
RE_DH = re.compile(r"(?:Delta\s*H|dH|deltaH)\s*[=:]\s*([-+0-9.eEdD]+)")
RE_PLAQ = re.compile(r"(?:[Aa]verage\s+plaquette|[Pp]laquette)\s*[=:]\s*([-+0-9.eEdD]+)")
RE_ACCEPT = re.compile(r"\b(Accept|Reject)\b", re.IGNORECASE)
RE_ACCRATE = re.compile(r"[Aa]cceptance\s*[=:]\s*([-+0-9.eEdD]+)")


def _f(s):
    return float(s.replace("d", "e").replace("D", "e"))


def parse_log(path):
    """Return rows as dict-of-arrays with keys traj, dH, acc, plaq (any may be empty).

    Strategy: a new 'Trajectory N' line flushes the current record. This matches the
    common Grid block layout. If almost no 'Trajectory' markers are found, fall back to
    anchoring one record per 'Delta H' line (one dH per MD trajectory).
    """
    rows, cur = [], {}
    n_traj_markers = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            mt = RE_TRAJ.search(line)
            if mt:
                n_traj_markers += 1
                if cur:
                    rows.append(cur)
                    cur = {}
                cur["traj"] = int(mt.group(1))
            md = RE_DH.search(line)
            if md:
                cur["dH"] = _f(md.group(1))
            mp = RE_PLAQ.search(line)
            if mp:
                cur["plaq"] = _f(mp.group(1))
            ma = RE_ACCEPT.search(line)
            if ma:
                cur["acc"] = 1.0 if ma.group(1).lower() == "accept" else 0.0
            mr = RE_ACCRATE.search(line)
            if mr:
                cur["accrate"] = _f(mr.group(1))
    if cur:
        rows.append(cur)

    # Fallback: no useful trajectory markers -> one row per dH occurrence.
    if n_traj_markers < 2:
        rows, cur = [], {}
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                md = RE_DH.search(line)
                if md:
                    if cur:
                        rows.append(cur)
                    cur = {"dH": _f(md.group(1))}
                    continue
                mp = RE_PLAQ.search(line)
                if mp:
                    cur["plaq"] = _f(mp.group(1))
                ma = RE_ACCEPT.search(line)
                if ma:
                    cur["acc"] = 1.0 if ma.group(1).lower() == "accept" else 0.0
        if cur:
            rows.append(cur)
    return rows


def parse_table(path):
    a = np.loadtxt(path, ndmin=2)
    keys = ["traj", "dH", "acc", "plaq"]
    return [dict(zip(keys, row)) for row in a[:, :4]]


def _col(rows, key):
    return np.array([r[key] for r in rows if key in r], dtype=float)


def jackknife(x, func=np.mean):
    """Jackknife mean and error of func over samples x (1D)."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 2:
        return (func(x) if n else np.nan), np.nan
    full = func(x)
    drop = np.array([func(np.delete(x, i)) for i in range(n)])
    err = np.sqrt((n - 1) / n * np.sum((drop - drop.mean()) ** 2))
    return full, err


def tau_int(series, c=5.0):
    """Madras-Sokal integrated autocorrelation time with automatic windowing."""
    x = np.asarray(series, float)
    n = len(x)
    if n < 4:
        return np.nan, np.nan
    x = x - x.mean()
    var = np.dot(x, x) / n
    if var == 0:
        return 0.5, np.nan
    rho = np.array([np.dot(x[: n - t], x[t:]) / (n * var) for t in range(n)])
    tau = 0.5
    W = 1
    for W in range(1, n):
        tau = 0.5 + np.sum(rho[1 : W + 1])
        if W >= c * tau:
            break
    err = tau * np.sqrt((2.0 * (2 * W + 1)) / n)  # Madras-Sokal error estimate
    return tau, err


def drift_test(x):
    """Compare first-half vs second-half means; return (mean1, mean2, n_sigma)."""
    x = np.asarray(x, float)
    h = len(x) // 2
    a, b = x[:h], x[h:]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    m1, e1 = a.mean(), a.std(ddof=1) / np.sqrt(len(a))
    m2, e2 = b.mean(), b.std(ddof=1) / np.sqrt(len(b))
    nsig = abs(m1 - m2) / np.hypot(e1, e2) if np.hypot(e1, e2) > 0 else np.nan
    return m1, m2, nsig


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("log", help="Grid HMC log file (or a table with --table)")
    p.add_argument("--table", action="store_true",
                   help="parse as a whitespace table: columns traj dH acc plaq")
    p.add_argument("--therm", type=int, default=None,
                   help="thermalisation cut: # of leading trajectories to drop "
                        "(default: auto-estimate from the plaquette running mean)")
    p.add_argument("--out", default="sims/output/lattice", help="output dir for plots")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    rows = parse_table(args.log) if args.table else parse_log(args.log)
    if not rows:
        print("ERROR: parsed 0 records. Paste me ~20 lines of the log and I'll fix the "
              "regexes (or use --table).", file=sys.stderr)
        return 1

    traj = _col(rows, "traj")
    dH = _col(rows, "dH")
    acc = _col(rows, "acc")
    plaq = _col(rows, "plaq")
    n = len(rows)
    if traj.size != n:
        traj = np.arange(n, dtype=float)

    print(f"\nParsed {n} trajectory records from {args.log}")
    print(f"  fields populated:  dH={dH.size}  acc={acc.size}  plaq={plaq.size}  "
          f"traj-id={_col(rows, 'traj').size}")

    warn = False

    # ---- thermalisation cut --------------------------------------------------------------
    if args.therm is not None:
        cut = args.therm
    elif plaq.size > 20:
        run = np.cumsum(plaq) / np.arange(1, plaq.size + 1)
        final = plaq[plaq.size // 2:].mean()
        band = 0.5 * plaq[plaq.size // 2:].std(ddof=1)
        inside = np.where(np.abs(run - final) <= max(band, 1e-9))[0]
        cut = int(inside[0]) if inside.size else plaq.size // 5
    else:
        cut = 0
    cut = max(0, min(cut, n - 2))
    print(f"  thermalisation cut: drop first {cut} trajectories "
          f"({'manual' if args.therm is not None else 'auto'})")

    def after(x):
        return x[cut:] if x.size == n else x[cut:] if x.size > cut else x

    # ---- 1) acceptance -------------------------------------------------------------------
    print("\n--- HMC health -------------------------------------------------")
    if acc.size:
        a = after(acc)
        rate = a.mean()
        flag = "OK  " if 0.65 <= rate <= 0.92 else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[{flag}] acceptance        = {rate:6.3f}   (want 0.70-0.90)")
    else:
        print("[--- ] acceptance        = (no Accept/Reject lines found)")

    # ---- 2) <exp(-dH)> == 1  (the unbiasedness check) ------------------------------------
    if dH.size:
        d = after(dH)
        em, ee = jackknife(np.exp(-d))
        nsig = abs(em - 1.0) / ee if ee and np.isfinite(ee) else np.nan
        flag = "OK  " if (np.isfinite(nsig) and nsig < 3) else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[{flag}] <exp(-dH)>        = {em:6.4f} +/- {ee:6.4f}   "
              f"({nsig:.1f}sigma from 1; want <3)")
        dm, de = jackknife(d)
        print(f"       <dH>              = {dm:6.4f} +/- {de:6.4f}   "
              f"(max |dH| = {np.abs(d).max():.2f})")
    else:
        print("[--- ] <exp(-dH)>        = (no Delta H lines found)")

    # ---- 3) plaquette: plateau, drift, autocorrelation -----------------------------------
    if plaq.size:
        pl = plaq[cut:] if plaq.size == n else plaq[cut:]
        pm, pe = jackknife(pl)
        m1, m2, nsig = drift_test(pl)
        tau, terr = tau_int(pl)
        flag = "OK  " if (np.isfinite(nsig) and nsig < 3) else "WARN"
        warn |= flag.strip() == "WARN"
        print(f"[{flag}] plaquette (therm) = {pm:8.6f} +/- {pe:8.6f}")
        print(f"       drift 1st vs 2nd half: {m1:.6f} vs {m2:.6f}  ({nsig:.1f}sigma; want <3)")
        if np.isfinite(tau):
            print(f"       tau_int(plaq)     = {tau:5.2f} +/- {terr:4.2f} traj  "
                  f"-> ~{max(1, int(round(2 * tau)))} traj between independent configs")
            print(f"       N_eff (post-therm) ~ {pl.size / (2 * tau):.0f} independent samples")
    else:
        print("[--- ] plaquette         = (no plaquette lines found)")

    print("----------------------------------------------------------------")
    print(("VERDICT: WARN -- check the flagged line(s) above."
           if warn else "VERDICT: all checks pass. Run looks healthy."))

    # ---- plots ---------------------------------------------------------------------------
    if not args.no_plots:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover
            print(f"(plots skipped: {e})")
            return 1 if warn else 0
        os.makedirs(args.out, exist_ok=True)
        npan = sum(x.size > 0 for x in (plaq, dH))
        fig, axes = plt.subplots(max(1, npan + (dH.size > 0)), 1,
                                 figsize=(9, 2.6 * max(1, npan + (dH.size > 0))))
        axes = np.atleast_1d(axes)
        k = 0
        if plaq.size:
            ax = axes[k]; k += 1
            ax.plot(traj[: plaq.size], plaq, lw=0.8, color="C0")
            ax.axvline(traj[cut] if cut < traj.size else cut, color="0.5", ls="--",
                       label=f"therm cut ({cut})")
            ax.set_ylabel("plaquette"); ax.legend(fontsize=8); ax.grid(alpha=0.2)
            ax.set_title("HMC health: plaquette history")
        if dH.size:
            ax = axes[k]; k += 1
            ax.plot(traj[: dH.size], dH, lw=0.7, color="C3")
            ax.axhline(0, color="0.6", lw=0.6)
            ax.set_ylabel("dH"); ax.grid(alpha=0.2)
            ax = axes[k]; k += 1
            ax.hist(after(dH), bins=40, color="C3", alpha=0.8)
            ax.set_xlabel("dH"); ax.set_ylabel("count")
            ax.set_title(f"dH distribution (post-therm, <exp(-dH)>={em:.3f})")
        axes[-1].set_xlabel("trajectory")
        fig.tight_layout()
        out = os.path.join(args.out, "hmc_health.png")
        fig.savefig(out, dpi=130); plt.close(fig)
        print(f"\nwrote {out}")

    return 1 if warn else 0


if __name__ == "__main__":
    raise SystemExit(main())
