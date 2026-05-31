#!/usr/bin/env python3
r"""Where does our universe sit in the supraverse generation tree?

Left: the generation-depth distribution P(n) of universes for sub-, near-, and
super-critical branching ratios m (mean viable children per universe). The shape
of the foam's depth distribution is controlled by m, not assumed.

Right: the observer-weighted posterior on OUR depth. The decisive datum -- we
observe dark matter AND dark energy, both of which require a parent in this
framework -- zeroes the n=0 (OGU) bin: we are categorically not an original-
generation universe. Given that floor, the branching prior says whether we are
shallow or deep.

Writes figures/pdf/generation_depth.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import population as pop

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    D = 12
    regimes = [(0.6, "C0", "subcritical $m=0.6$ (lineage dies out)"),
               (1.0, "C1", "critical $m=1$ (scale-free)"),
               (1.8, "C3", "supercritical $m=1.8$ (deep tree)")]

    lines = ["Supraverse generation-depth distribution",
             "=" * 42,
             "Branching ratio m = mean VIABLE children per universe",
             "(super-critical AND past the chiral-vortical birth filter).",
             ""]
    for m, _, lab in regimes:
        en = pop.expected_generation(m, D)
        lines.append(f"  m={m:<4} <n> (universe-counted) = {en:5.2f}")
    lines += [
        "",
        "ROBUST ANCHOR (independent of m):",
        "  We observe dark matter AND dark energy.",
        "  Both require a parent (DM = projected parent matter; DE = parent",
        "  accretion). An OGU (n=0) has no parent -> would have neither.",
        "  => Our universe is NOT an OGU.  n >= 1.",
        "",
        f"  we_are_ogu(False,False) = {pop.we_are_ogu(False, False)}",
        f"  we_are_ogu(True,True)   = {pop.we_are_ogu(True, True)}",
        f"  min generation from observations = "
        f"{pop.min_generation_from_observations(True, True)}",
    ]
    # posterior depth for the supercritical case, observers slightly decaying
    post = pop.depth_posterior(m=1.8, D=D, n_min=1, structure_decay=0.7)
    n = np.arange(D + 1)
    n_map = int(n[np.argmax(post)])
    lines += ["",
              f"Illustrative posterior (m=1.8, structure_decay=0.7, n>=1):",
              f"  most probable depth n* = {n_map}, "
              f"<n> = {(n*post).sum():.2f}",
              "  (parameters illustrative; tightened by DM/DE amplitudes + C.)"]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "generation_depth.txt"), "w") as f:
        f.write(text + "\n")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # --- Panel 1: depth distribution vs branching ratio ---
    for m, c, lab in regimes:
        p = pop.generation_pmf(m, D)
        ax1.plot(n, p, "o-", color=c, ms=4, lw=1.4, label=lab)
    ax1.set_xlabel("generation depth $n$  (0 = OGU)")
    ax1.set_ylabel("fraction of universes $P(n)$")
    ax1.set_title("Foam depth distribution is set by branching ratio $m$")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.2)

    # --- Panel 2: our posterior with the n>=1 floor ---
    prior = pop.observer_weighted_pmf(1.8, D, structure_decay=0.7)
    ax2.bar(n, prior, width=0.8, color="0.8", label="branching prior (obs-wtd)")
    ax2.bar(n, post, width=0.55, color="C2",
            label="posterior given $n\\geq1$ (DM+DE)")
    ax2.axvspan(-0.5, 0.5, color="C3", alpha=0.12)
    ax2.annotate("OGU ruled out\n(no parent $\\Rightarrow$ no DM/DE)",
                 xy=(0, prior[0]), xytext=(1.4, max(prior) * 0.78),
                 fontsize=8, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    ax2.set_xlabel("our generation depth $n$")
    ax2.set_ylabel("probability")
    ax2.set_title("We are not an OGU: $n\\geq1$, depth a range for now")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.2)

    fig.suptitle("Birth vs growth: the supraverse generation tree and our place "
                 "in it", fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(PDF_DIR, "generation_depth.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
