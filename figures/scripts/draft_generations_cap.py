#!/usr/bin/env python3
r"""Generations from the S budget: the cap exists, requires walking, and 3 is sensitive.

The number of generations is capped because each adds to the electroweak S and the
total must stay under the precision budget (~0.1): N_max = budget / S_per_gen.

The plot: N_max vs S_per_generation (the steep 1/x cap). Marked:
 - leading-order S ~ 0.16 -> cap 0 (would forbid even our existence -> walking is
   MANDATORY, not optional);
 - the narrow window S_per_gen in [0.025, 0.033] that yields exactly 3;
 - the ~5x walk-down our three generations demand -- which the torsion's G_A=G_V
   (Fierz + RPA) supplies in the right direction.

Renders figures/pdf/generations_cap.pdf, writes sims/output/generations_cap.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import generations_cap as gc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    w = gc.requires_walking()
    lo3, hi3 = gc.s_per_gen_for_cap(3)
    lines = [
        "Generations from the electroweak S budget (cap = budget / S_per_gen)",
        "=" * 66,
        f"  leading-order S_per_gen = N_c/6pi = {gc.S_LEADING:.3f}  ->  cap = "
        f"{gc.n_max(gc.S_LEADING)}",
        "    (cap 0 would forbid even our own existence -> WALKING IS MANDATORY)",
        f"  to allow our 3 generations: S_per_gen <= {gc.S_BUDGET/3:.3f}  "
        f"(a {w['walk_factor_needed']:.0f}x walk-down from leading)",
        f"  cap = 3 exactly: S_per_gen in ({lo3:.3f}, {hi3:.3f}]",
        "",
        "  the cap is a steep 1/x law (highly sensitive):",
        "    S=0.05 -> 2,  S=0.033 -> 3,  S=0.025 -> 4,  S=0.016 -> 6",
        "",
        "READING: a cap EXISTS (the SM does not explain why generations are finite; here",
        "it follows from the S budget). Our existence DEMANDS walking (leading order",
        "forbids us). The torsion's G_A=G_V supplies walking in the right direction",
        "(Fierz + RPA). Three is squarely in the plausible window. But N_max is",
        "exponentially sensitive to S_per_gen, which is not yet pinned -- so the cap is",
        "real and three is plausible, but the exact count is owed (full chiral RPA).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "generations_cap.txt"), "w") as f:
        f.write(text + "\n")

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    s = np.logspace(np.log10(0.008), np.log10(0.25), 500)
    ax.step(s, [gc.n_max(x) for x in s], where="mid", color="C0", lw=2.0)

    # leading order: cap 0, forbids existence
    ax.axvline(gc.S_LEADING, color="r", lw=1.5, ls="--")
    ax.annotate("leading order\n$N_c/6\\pi\\approx0.16$\ncap = 0 (forbids us!)",
                (gc.S_LEADING, 4.5), fontsize=8.5, color="r", ha="center")
    # cap-3 window
    ax.axvspan(lo3, hi3, color="C2", alpha=0.25)
    ax.annotate("cap = 3 window\n$S\\in[0.025,0.033]$", (np.sqrt(lo3 * hi3), 7.5),
                fontsize=8.5, color="C2", ha="center")
    # the walk arrow
    ax.annotate("", xy=(0.03, 0.5), xytext=(0.15, 0.5),
                arrowprops=dict(arrowstyle="->", color="C1", lw=2.0))
    ax.text(0.062, 0.9, "torsion walking\n(Fierz: $G_A=G_V$)\nmust carry us here",
            fontsize=8.5, color="C1", ha="center")

    ax.set_xscale("log")
    ax.set_xlim(0.008, 0.25)
    ax.set_ylim(0, 11)
    ax.set_xlabel(r"$S$ per generation (walked; not yet pinned)")
    ax.set_ylabel(r"max number of generations  $N_{\max}=\lfloor 0.1/S\rfloor$")
    ax.set_title("Generations are capped by the S budget\n"
                 "the cap exists and requires walking; '3' is plausible but sensitive",
                 fontsize=12)
    ax.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "generations_cap.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "generations_cap.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'generations_cap.pdf')}")


if __name__ == "__main__":
    main()
