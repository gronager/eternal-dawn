#!/usr/bin/env python3
r"""Tier 1b: Oppenheimer--Snyder collapse with a torsion-bounce interior.

One figure, the whole thesis: from outside, the surface freezes at the horizon
and reddens to black (a black hole); from inside, the dust bounces at rho_C deep
behind the horizon and re-expands (an inverse bubble). Renders
figures/pdf/os_collapse.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import os_collapse as osc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    s = osc.simulate_os(a_init=8.0, rs_units=3.0)
    r = s.realistic

    lines = [
        "Tier 1b: Oppenheimer-Snyder collapse with a torsion bounce",
        "=" * 60,
        "INTERIOR (dust, w=0): bounces at R_min instead of a=0 singularity.",
        "EXTERIOR (Schwarzschild): surface freezes at the horizon, light -> black.",
        "",
        f"  realistic mass M           = {r['M_kg']:.3e} kg (observable universe)",
        f"  Schwarzschild radius r_s   = {r['r_s_m']:.3e} m",
        f"  bounce areal radius R_min  = {r['r_min_m']:.3f} m",
        f"  bounce proper volume       = {r['bounce_volume_m3']:.1f} m^3",
        f"  R_min / r_s (universe)     = {r['rmin_over_rs']:.2e}",
        f"  R_min / r_s (10 Msun BH)   = {r['rmin_over_rs_stellar']:.2e}",
        "",
        "  => the bounce is hidden FAR inside the horizon (Tardis effect):",
        "     ~6 m across at the bounce, ~10^26 m horizon, ~46 Gly once expanded.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "os_collapse.txt"), "w") as f:
        f.write(text + "\n")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))

    # Panel A: interior areal radius (log), bounce vs singularity.
    ax1.semilogy(s.tau, s.R_int, color="C0", lw=2.0,
                 label="Einstein--Cartan interior (bounce)")
    pos = s.R_os > 0
    ax1.semilogy(s.tau_os[pos], s.R_os[pos], color="C3", lw=1.6, ls="--",
                 label="general relativity (singularity)")
    ax1.axhline(s.r_s_units, color="0.4", lw=1.2)
    ax1.text(s.tau[-1], s.r_s_units * 1.15, "event horizon $r_s$",
             ha="right", color="0.3", fontsize=9)
    i_b = int(np.argmin(s.R_int))
    ax1.plot([s.tau[i_b]], [s.R_int[i_b]], "o", color="C0", ms=6)
    ax1.annotate(r"bounce $R_{\min}$ ($\rho=\rho_C$)", (s.tau[i_b], s.R_int[i_b]),
                 textcoords="offset points", xytext=(8, 6), fontsize=9)
    ax1.annotate(r"$a\to0$", (s.tau_os[pos][-1], s.R_os[pos][-1]),
                 textcoords="offset points", xytext=(0, -14), color="C3",
                 fontsize=9, ha="center")
    ax1.set_xlabel(r"interior proper time  $\tau$")
    ax1.set_ylabel(r"areal radius  $R / R_{\min}$")
    ax1.set_title("Interior: bounce hidden behind the horizon")
    ax1.legend(fontsize=8, loc="lower right")
    ax1.grid(True, which="both", alpha=0.2)

    # Panel B: exterior view -- frozen star.
    ax2.plot(s.R_ext, s.t_ext, color="C4", lw=2.0, label="Schwarzschild time $t$")
    ax2.axvline(1.0, color="0.4", lw=1.2)
    ax2.text(1.02, s.t_ext.max() * 0.6, "horizon", color="0.3", fontsize=9)
    ax2.set_xlabel(r"areal radius  $R / r_s$")
    ax2.set_ylabel(r"exterior time  $t / (r_s/c)$", color="C4")
    ax2.tick_params(axis="y", labelcolor="C4")
    ax2.set_title("Exterior: surface freezes, light reddens to black")
    ax2.invert_xaxis()                     # infall goes left, toward the horizon
    ax2b = ax2.twinx()
    ax2b.plot(s.R_ext, s.z_ext, color="C1", lw=1.6, ls="-.",
              label="redshift $1+z$")
    ax2b.set_ylabel(r"surface redshift  $1+z$", color="C1")
    ax2b.tick_params(axis="y", labelcolor="C1")
    ax2b.set_ylim(0, None)
    lines1, lab1 = ax2.get_legend_handles_labels()
    lines2, lab2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, lab1 + lab2, fontsize=8, loc="upper left")

    fig.suptitle("Oppenheimer--Snyder collapse with a torsion bounce: "
                 "black hole outside, bouncing cosmology inside", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_pdf = os.path.join(PDF_DIR, "os_collapse.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
