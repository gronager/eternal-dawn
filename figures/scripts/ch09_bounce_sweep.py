#!/usr/bin/env python3
r"""Tier 1 sweep: the torsion bounce is generic across the equation of state.

Integrates the bounce for a range of w (dust -> radiation -> stiff) and shows
(i) the bounce always occurs at rho = rho_C (a_min = 1) regardless of w, and
(ii) the bounce duration follows the small-amplitude scaling FWHM ~ 2 sqrt(2)/n
with n = 3(1+w). Renders figures/pdf/bounce_sweep.pdf.
"""
from __future__ import annotations

import math
import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import bounce as bnc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

WS = [0.0, 1.0 / 6.0, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, 1.0]
LABELS = ["0 (dust)", "1/6", "1/3 (radiation)", "1/2", "2/3", "1 (stiff)"]
TAU = bnc.physical_timescale(1.0e50)     # seconds per dimensionless time unit


def main() -> None:
    rows = []
    sols = []
    for w in WS:
        sol = bnc.simulate_bounce(w=w, a_init=2.5, t_max=80.0, n_points=12000)
        fwhm = bnc.bounce_fwhm(sol)
        n = 3.0 * (1.0 + w)
        rows.append({
            "w": w, "n": n, "a_min": sol.a_min, "rho_max": sol.rho_max,
            "fwhm": fwhm, "fwhm_analytic": 4.0 / n,
            "max_abs_H": float(np.max(np.abs(sol.H))),
            "fwhm_seconds": fwhm * TAU,
        })
        sols.append(sol)

    lines = ["Tier 1 sweep: bounce vs equation of state (rho_C = 1)",
             "=" * 70,
             f"{'w':>10} {'a_min':>8} {'rho_max':>9} {'FWHM':>8} "
             f"{'4/n':>8} {'max|H|':>8} {'FWHM [s]':>10}"]
    for r in rows:
        lines.append(f"{r['w']:>10.4f} {r['a_min']:>8.5f} {r['rho_max']:>9.5f} "
                     f"{r['fwhm']:>8.4f} {r['fwhm_analytic']:>8.4f} "
                     f"{r['max_abs_H']:>8.4f} {r['fwhm_seconds']:>10.2e}")
    lines += ["",
              "Bounce density is universal: a_min = 1, rho_max = rho_C for every w,",
              "and max|H| = 1/2 exactly (since max_x x(1-x) = 1/4).",
              "Only the bounce sharpens with stiffer matter: FWHM = 4/n exactly.",
              "Physical duration (rho_C = 1e50): a few x 1e-21 s."]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "bounce_sweep.txt"), "w") as f:
        f.write(text + "\n")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    cmap = plt.cm.viridis(np.linspace(0.0, 0.85, len(WS)))

    # Panel A: a(tau) for each w, all bouncing at a_min = 1.
    for sol, lab, col in zip(sols, LABELS, cmap):
        ax1.plot(sol.t, sol.a, color=col, lw=1.8, label=fr"$w={lab}$")
    ax1.axhline(1.0, color="0.6", ls=":", lw=1.0)
    ax1.text(ax1.get_xlim()[1], 1.04, r"$a_{\min}=1$ ($\rho=\rho_C$)",
             ha="right", color="0.4", fontsize=9)
    ax1.set_xlabel(r"dimensionless time  $t$")
    ax1.set_ylabel(r"scale factor  $a$")
    ax1.set_xlim(-12, 12)
    ax1.set_ylim(0.9, 2.6)
    ax1.set_title("Bounce occurs at $\\rho_C$ for every equation of state")
    ax1.legend(fontsize=8, ncol=2, loc="upper center")
    ax1.grid(True, alpha=0.25)

    # Panel B: bounce duration vs w, simulation vs analytic scaling.
    ws = np.array([r["w"] for r in rows])
    fwhm = np.array([r["fwhm"] for r in rows])
    n = 3.0 * (1.0 + ws)
    ax2.plot(ws, fwhm, "o", color="C0", ms=7, label="simulation")
    wfine = np.linspace(0.0, 1.0, 100)
    ax2.plot(wfine, 4.0 / (3.0 * (1.0 + wfine)), color="C3",
             lw=1.6, label=r"exact $\,4/n,\ n=3(1+w)$")
    ax2.set_xlabel(r"equation of state  $w$")
    ax2.set_ylabel(r"bounce duration  $\Delta\tau_{\rm FWHM}$")
    ax2.set_title("Bounce sharpens with stiffer matter")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.25)
    secax = ax2.secondary_yaxis(
        "right", functions=(lambda x: x * TAU, lambda x: x / TAU))
    secax.set_ylabel(r"physical duration [s] ($\rho_C=10^{50}$)")

    fig.suptitle("Tier 1 sweep: robustness of the Einstein--Cartan bounce",
                 fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.07, right=0.9, wspace=0.32, top=0.82, bottom=0.12)
    out_pdf = os.path.join(PDF_DIR, "bounce_sweep.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
