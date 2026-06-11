#!/usr/bin/env python3
r"""Paper C, Fig. 2 -- the pion and torsiton (nucleon) effective-mass plateaus.

The binding result of Sec. III.B: a clean effective-mass plateau at m_N > m_pi > 0 is the torsiton,
found non-perturbatively. m_eff(t) = ln[C(t)/C(t+1)] for the pion (calibration) and the nucleon (the
torsiton), with delete-1 jackknife errors over configs and the plateau-fit bands. Renders
figures/pdf/torsiton_plateaus.pdf.

Data: figures/data/baryon_L16x48.dat -- the run/06 output (rows: cfg t C_pi C_N) copied from
out/dyn_L16x48_m-0.5/baryon_raw.dat. (The overall sign of C_N is a contraction convention and cancels
in the ratio.) If the file is absent the script says so and exits cleanly.
"""
from __future__ import annotations

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import lattice as lat

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
DATA = os.path.join(ROOT, "figures", "data", "baryon_L16x48.dat")
T = 48
os.makedirs(PDF_DIR, exist_ok=True)


def _jackknife_meff(C):
    """C: (ncfg, nt) correlator per config. Returns t, meff(t), meff_err(t) (delete-1 jackknife).
    m_eff(t) = ln[Cbar(t)/Cbar(t+1)]; the ratio is sign-robust (an overall sign cancels)."""
    ncfg, nt = C.shape
    full = C.mean(0)
    meff_full = np.log(full[:-1] / full[1:])
    jk = np.empty((ncfg, nt - 1))
    for i in range(ncfg):
        m = np.delete(C, i, axis=0).mean(0)
        jk[i] = np.log(m[:-1] / m[1:])
    err = np.sqrt((ncfg - 1) / ncfg * np.sum((jk - jk.mean(0)) ** 2, axis=0))
    return np.arange(nt - 1), meff_full, err


def main():
    if not os.path.exists(DATA):
        print(f"missing {DATA}\n  -> copy your run/06 output in, e.g.:\n"
              f"     cp out/dyn_L16x48_m-0.5/baryon_raw.dat figures/data/baryon_L16x48.dat\n"
              f"  then re-run this script.")
        sys.exit(0)

    raw = np.loadtxt(DATA)
    if raw.ndim != 2 or raw.shape[1] != 4:
        sys.exit("expected rows 'cfg t C_pi C_N'")
    cfgs = np.unique(raw[:, 0]).astype(int)
    ts = np.unique(raw[:, 1]).astype(int)
    nt = ts.max() + 1
    Cpi = np.full((len(cfgs), nt), np.nan)
    CN = np.full((len(cfgs), nt), np.nan)
    for i, c in enumerate(cfgs):
        sub = raw[raw[:, 0] == c]
        for t, cp, cn in sub[:, 1:]:
            Cpi[i, int(t)] = cp
            CN[i, int(t)] = cn

    tpi, mpi, mpi_e = _jackknife_meff(Cpi)
    tN, mN, mN_e = _jackknife_meff(CN)

    # the published plateau fits (from cartasis_sims.lattice.baryon_spectrum on this ensemble)
    res = lat.baryon_spectrum(raw, T=T)
    fits = {"pion": ("C0", res["pion"]), "nucleon": ("C3", res["nucleon"])}

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    for (tt, mm, ee, col, lab) in [
        (tpi, mpi, mpi_e, "C0", r"pion  $m_\pi$"),
        (tN, mN, mN_e, "C3", r"torsiton  $m_N$"),
    ]:
        ax.errorbar(tt, mm, yerr=ee, fmt="o", color=col, ms=4, lw=1.2, capsize=2.5, label=lab)
    for name, (col, r) in fits.items():
        if np.isfinite(r["mass"]):
            t0, t1, m, e = r["tmin"], r["tmax"], r["mass"], r["mass_err"]
            ax.fill_between([t0, t1], m - e, m + e, color=col, alpha=0.25)
            ax.annotate(fr"$m\,a={m:.3f}({int(round(e*1000))})$", (t1, m), xytext=(6, 0),
                        textcoords="offset points", color=col, fontsize=9, va="center")

    ax.set_xlim(0, T // 2)
    ax.set_ylim(0, None)
    ax.set_xlabel(r"$t/a$")
    ax.set_ylabel(r"effective mass  $m_{\rm eff}(t)\,a = \ln[\,C(t)/C(t{+}1)\,]$")
    ax.set_title("Pion and torsiton effective-mass plateaus\n"
                 r"$m_N > m_\pi > 0$: the torsiton binds ($m_N/m_\pi\to 3/2$, constituent counting)",
                 fontsize=11)
    ax.legend(fontsize=9.5, loc="upper right")
    ax.grid(alpha=0.18)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "torsiton_plateaus.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")
    if np.isfinite(res["pion"]["mass"]) and np.isfinite(res["nucleon"]["mass"]):
        print(f"  m_pi a = {res['pion']['mass']:.3f}, m_N a = {res['nucleon']['mass']:.3f}, "
              f"m_N/m_pi = {res['nucleon']['mass']/res['pion']['mass']:.3f}")


if __name__ == "__main__":
    main()
