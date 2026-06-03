#!/usr/bin/env python3
r"""The bounce as a potential wall (Ch.2).

Writing the effective Friedmann equation H^2 = a^{-n}(1 - a^{-n}) as a zero-energy
1-D mechanics problem, a'^2 + V(a) = 0, gives the effective potential

    V(a) = a^{2-2n} - a^{2-n},   n = 3(1+w).

The TORSION term (a^{2-2n}, from the -rho^2/rho_C spin contribution) is a repulsive
wall at small a; the MATTER term (-a^{2-n}) is the attractive pull at large a; they
cross at a = a_min, where the contracting universe turns around. It is a SINGLE
repulsive wall -- not a multi-well potential -- whose only knobs are the equation of
state w and the density scale rho_C. Renders figures/pdf/bounce_potential.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import bounce as bn

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    a = np.linspace(0.78, 2.6, 600)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))

    for w, lab, c in [(0.0, "matter  $w=0$  ($n=3$)", "C0"),
                      (1 / 3, "radiation  $w=1/3$  ($n=4$)", "C1"),
                      (1.0, "stiff  $w=1$  ($n=6$)", "C2")]:
        V = bn.effective_potential(a, w)
        ax.plot(a, V, color=c, lw=2.0, label=lab)
        # decompose for the radiation case to show the two competing terms
        if abs(w - 1 / 3) < 1e-9:
            n = 3 * (1 + w)
            ax.plot(a, a**(2 - 2 * n), color=c, lw=1.0, ls=":",
                    alpha=0.7, label="  torsion wall  $a^{2-2n}$")
            ax.plot(a, -a**(2 - n), color=c, lw=1.0, ls="--",
                    alpha=0.7, label="  matter pull  $-a^{2-n}$")

    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(1.0, color="0.5", lw=1.0, ls=":")
    ax.plot(1.0, 0.0, "ko", ms=7)
    ax.annotate("bounce\n$a=a_{\\min}$, $V=0$", (1.0, 0.0),
                textcoords="offset points", xytext=(10, 14), fontsize=9)
    ax.annotate("forbidden\n($V>0$)", (0.85, 0.5), fontsize=8.5, color="0.4")
    ax.annotate("allowed ($V<0$): the expanding universe", (1.5, -0.55),
                fontsize=8.5, color="0.4")

    ax.set_xlim(0.78, 2.6)
    ax.set_ylim(-0.9, 1.2)
    ax.set_xlabel(r"scale factor  $a/a_{\min}$")
    ax.set_ylabel(r"effective potential  $V(a)$  ($\dot a^2 + V = 0$)")
    ax.set_title("The bounce is a single repulsive wall\n"
                 "torsion ($a^{2-2n}$) vs matter ($-a^{2-n}$); the only knobs are "
                 "$w$ and $\\rho_C$", fontsize=11.5)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "bounce_potential.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "bounce_potential.png"), dpi=130)
    plt.close(fig)
    print(f"wrote {os.path.join(PDF_DIR, 'bounce_potential.pdf')}")


if __name__ == "__main__":
    main()
