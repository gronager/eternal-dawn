#!/usr/bin/env python3
r"""Tier 1 (minimal): the torsion bounce vs the GR singularity.

Integrates a collapsing radiation-dominated spin-fluid universe through the
Einstein--Cartan bounce and contrasts it with the same collapse in plain GR,
which reaches a = 0 (a singularity) in finite time. Renders
figures/pdf/bounce.pdf.
"""
from __future__ import annotations

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


def main() -> None:
    w = 1.0 / 3.0          # radiation: matter is ultra-relativistic at rho_C
    a_init = 3.0
    sol = bnc.simulate_bounce(w=w, a_init=a_init, t_max=14.0)
    gr = bnc.simulate_gr_collapse(w=w, a_init=a_init)

    # Align the GR collapse so it starts at the same (t, a) as the EC collapse.
    t_start = sol.t[0]                       # negative: start of collapse branch
    gr_t = gr["t"] + t_start
    tau = bnc.physical_timescale(1.0e50)     # seconds per dimensionless unit

    summary = [
        "Tier 1: Einstein-Cartan torsion bounce",
        "=" * 46,
        f"  equation of state w        = {w:.3f} (radiation)",
        f"  scale factor at bounce a_min = {sol.a_min:.6f}  (target 1.0)",
        f"  max density rho_max/rho_C  = {sol.rho_max:.6f}  (target 1.0)",
        f"  H at bounce                = {sol.H[np.argmin(sol.a)]:.2e}  (target 0)",
        f"  GR collapse reaches a=0 at t = {gr['t_singularity']:.3f} (singularity)",
        f"  physical time unit (rho_C=1e50) = {tau:.2e} s",
        f"  => bounce duration ~ {tau:.0e} s",
    ]
    text = "\n".join(summary)
    print(text)
    with open(os.path.join(OUT_DIR, "bounce.txt"), "w") as f:
        f.write(text + "\n")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: scale factor, EC bounce vs GR singularity.
    ax1.plot(sol.t, sol.a, color="C0", lw=2.0, label="Einstein--Cartan (bounce)")
    ax1.plot(gr_t, gr["a"], color="C3", lw=1.6, ls="--",
             label="general relativity (singularity)")
    ax1.axhline(sol.a_min, color="0.6", ls=":", lw=1.0)
    ax1.plot([0], [sol.a_min], "o", color="C0", ms=6)
    ax1.annotate("bounce at $a_{\\min}$ ($\\rho=\\rho_C$)", (0, sol.a_min),
                 textcoords="offset points", xytext=(10, 6), fontsize=9)
    ax1.plot([gr_t[-1]], [gr["a"][-1]], "x", color="C3", ms=8)
    ax1.annotate("$a\\to0$", (gr_t[-1], gr["a"][-1]),
                 textcoords="offset points", xytext=(-30, 8), color="C3",
                 fontsize=9)
    ax1.set_xlabel("dimensionless time  $t$")
    ax1.set_ylabel("scale factor  $a$")
    ax1.set_ylim(0, a_init * 1.05)
    ax1.set_title("Collapse: bounce vs singularity")
    ax1.legend(fontsize=8, loc="upper center")
    ax1.grid(True, alpha=0.25)

    # Panel B: density and Hubble rate through the bounce.
    ax2.plot(sol.t, sol.rho, color="C1", lw=2.0, label=r"$\rho/\rho_C$")
    ax2.axhline(1.0, color="0.6", ls=":", lw=1.0)
    ax2.set_xlabel("dimensionless time  $t$")
    ax2.set_ylabel(r"density  $\rho/\rho_C$", color="C1")
    ax2.tick_params(axis="y", labelcolor="C1")
    ax2.set_ylim(0, 1.15)
    ax2b = ax2.twinx()
    ax2b.plot(sol.t, sol.H, color="C2", lw=1.6, label="$H$")
    ax2b.axhline(0.0, color="C2", ls=":", lw=0.8)
    ax2b.set_ylabel("Hubble rate  $H$", color="C2")
    ax2b.tick_params(axis="y", labelcolor="C2")
    ax2.set_title(r"Density saturates at $\rho_C$, $H$ crosses zero")
    ax2.set_xlim(sol.t[0], sol.t[-1])

    fig.suptitle("Einstein--Cartan torsion bounce (Tier 1, radiation fluid)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_pdf = os.path.join(PDF_DIR, "bounce.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
