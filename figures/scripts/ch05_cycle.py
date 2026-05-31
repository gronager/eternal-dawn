#!/usr/bin/env python3
r"""The recursive cycle: why a heat-dead universe can seed the next foam.

The total entropy rises monotonically (Second Law), but the Weyl (clumping)
gravitational entropy CYCLES -- low at the bounce, high while black holes exist, low
again once they evaporate and the interior smooths to de Sitter. A heat-dead
interior is thus a low-gravitational-entropy void, ready to nucleate the next foam:
Penrose's CCC, made to branch. Renders figures/pdf/cycle.pdf, writes
sims/output/cycle.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import cycle as cy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    H_lambda = 1.8e-18
    T_dS = cy.de_sitter_temperature(H_lambda)

    lines = [
        "The recursive cycle: a heat-dead universe is the next foam's void",
        "=" * 62,
        "  Penrose's CCC made to branch: the smooth, empty de Sitter far future is a",
        "  low-Weyl-entropy void, and it nucleates a whole foam of fresh OGUs.",
        "",
        "  The apparent paradox (heat death = max entropy, nucleation needs low",
        "  entropy) is resolved by Penrose's distinction:",
        f"    Weyl (clumping) entropy:  low -> high (black holes) -> low (de Sitter)",
        f"      bounce  ~10^{cy.weyl_entropy(cy.T_BOUNCE):.0f},  peak ~10^"
        f"{cy.weyl_entropy(3e40):.0f},  heat death ~10^{cy.weyl_entropy(cy.T_DESITTER):.0f}",
        f"    Total entropy:            monotonic  10^{cy.total_entropy(cy.T_BOUNCE):.0f}"
        f" -> 10^{cy.total_entropy(cy.T_DESITTER):.0f}",
        "",
        f"  de Sitter void temperature T_dS = {T_dS:.1e} K  (= the foam-coarsening",
        f"  bath, and ~T_H of a horizon-mass hole -- one temperature for the whole",
        f"  far future). The void lasts forever, so even a tiny nucleation rate fires:",
        f"  the cycle always closes.",
        "",
        "READING: yes -- it cycles. A ginormous heat-dead OGU is, gravitationally, a",
        "smooth empty void at ~1e-30 K, indistinguishable from the substrate SCT",
        "nucleates OGUs from -- because it IS that substrate. The void is not a",
        "primordial 'outside'; it is the asymptotic state every universe reaches. So",
        "universes nest both ways: we are inside a parent, and our own far future will",
        "breed a foam of fresh OGUs in its emptiness. No first cause, no last; the",
        "supraverse is self-similar and eternal, with the black-hole Hawking clock",
        "(~1e100 yr) setting the reset between aeons. This is the foam thread's",
        "natural terminus -- and Penrose's CCC with branching instead of a single heir.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "cycle.txt"), "w") as f:
        f.write(text + "\n")

    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    t = np.logspace(-36, 117, 2000)
    ax.plot(t, cy.total_entropy(t), "C0", lw=2.4,
            label="total entropy (thermal + horizon) -- monotonic (Second Law)")
    ax.plot(t, cy.weyl_entropy(t), "C3", lw=2.4,
            label="Weyl / clumping entropy -- low $\\to$ high $\\to$ low (resets)")

    # phase shading and labels
    ax.axvspan(1e-36, cy.T_STRUCTURE, color="C2", alpha=0.07)
    ax.axvspan(cy.T_STRUCTURE, cy.T_BH_EVAP, color="C1", alpha=0.07)
    ax.axvspan(cy.T_BH_EVAP, 1e117, color="C0", alpha=0.07)
    ax.text(1e-20, 70, "bounce\n(smooth,\nlow Weyl)", fontsize=8, color="C2",
            ha="center")
    ax.text(1e60, 108, "structure + black holes\n(high Weyl)", fontsize=8,
            color="C1", ha="center")
    ax.text(1e113, 70, "de Sitter heat death\n(smooth void,\nlow Weyl again)",
            fontsize=8, color="C0", ha="center")

    # the recursion arrow: heat-dead void -> nucleate the next bounce
    ax.annotate("", xy=(1e-30, 80), xytext=(2e115, 80),
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.6,
                                connectionstyle="arc3,rad=-0.35"))
    ax.text(1e40, 64, "the smooth void nucleates the next foam of OGUs "
            "(Penrose CCC, branching)", fontsize=8.5, color="0.35", ha="center")

    ax.set_xscale("log")
    ax.set_xlabel(r"time since bounce  [s]   ($\sim10^{108}$ s $=10^{100}$ yr: black-hole evaporation)")
    ax.set_ylabel(r"entropy  $\log_{10} S$  [$k_B$]")
    ax.set_title("The recursive cycle: total entropy rises, but Weyl entropy resets\n"
                 "so a heat-dead interior is a low-gravitational-entropy void for the "
                 "next foam", fontsize=11)
    ax.set_ylim(60, 126)
    ax.legend(fontsize=9, loc="center left")
    ax.grid(True, which="both", alpha=0.15)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "cycle.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
