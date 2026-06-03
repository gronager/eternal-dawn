#!/usr/bin/env python3
r"""Generations as an overlap hierarchy: arithmetic levels -> geometric masses.

A: the real Standard-Model charged-fermion masses, log scale, vs generation. They
   are approximately STRAIGHT lines -- evenly spaced in log, i.e. geometric -- the
   structure an exponential-overlap mechanism produces. A tuned soliton overlap
   hierarchy is overlaid.
B: the mechanism -- the overlap O_n of the n-th internal level with a localized core
   source falls EXPONENTIALLY with n (straight on a log axis), so mass ~ |O_n|^2 is
   geometric. Steeper source -> bigger per-generation ratio.

Renders figures/pdf/generations.pdf, writes sims/output/generations.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import generations as gen
from cartasis_sims import soliton as so

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    lv = so.energy_levels(n_levels=4, kind="linear", depth=6.0, E_max=18.0, n_scan=110)

    lines = ["Generations: arithmetic levels -> geometric masses (overlap hierarchy)",
             "=" * 68,
             "  Standard-Model fermion masses are approximately GEOMETRIC",
             "  (ln-mass gaps equal to within a factor ~2):"]
    for name, m in gen.SM_FERMIONS.items():
        g = gen.log_spacings(m)
        lines.append(f"    {name:22s} ln-gaps = {np.round(g,2)}  "
                     f"(ratio {max(g)/min(g):.1f}; geometric-ish={gen.is_approximately_geometric(m)})")
    lines += ["",
              "  the overlap mechanism (soliton levels, core source) -- tunable steepness:"]
    for rs in (0.15, 0.25, 0.4):
        m = gen.overlap_masses(lv, source_size=rs)
        lines.append(f"    source_size={rs}: per-generation factor = {gen.geometric_factor(m):.1f}x")
    lines += ["",
              "READING: an arithmetic ladder of internal levels, mapped through an",
              "EXPONENTIALLY-sensitive overlap, gives ln(mass) linear in the level index --",
              "a geometric hierarchy. The SM's 'unnaturally tiny' electron Yukawa (~3e-6) is",
              "then just exp(-O(10)), the exponential of an order-one overlap, not a tuned",
              "small number. The mechanism reproduces the SHAPE (geometric, orders of",
              "magnitude from O(1) inputs); the exact ratios (leptons 207, 17 -- a DECREASING",
              "ratio) need the detailed level/overlap structure and are the owed refinement."]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "generations.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.2))

    # ---- A: SM masses (log) vs generation, ~geometric ----
    gens = np.array([0, 1, 2])
    colors = {"leptons (e, mu, tau)": "C0", "up-type (u, c, t)": "C1",
              "down-type (d, s, b)": "C3"}
    for name, m in gen.SM_FERMIONS.items():
        axA.semilogy(gens, m, "o-", color=colors[name], lw=1.6, ms=7, label=name)
        # a straight (geometric) reference line through the 1st and 3rd of each family
        ref = m[0] * (m[2] / m[0]) ** (gens / 2.0)
        axA.semilogy(gens, ref, "--", color=colors[name], lw=0.8, alpha=0.5)
    axA.set_xticks(gens); axA.set_xticklabels(["1st", "2nd", "3rd"])
    axA.set_xlabel("generation")
    axA.set_ylabel("mass [MeV]")
    axA.set_title("SM fermion masses are ~geometric\n(evenly spaced in log; dashed = exact geometric)",
                  fontsize=11)
    axA.legend(fontsize=8.5, loc="lower right")
    axA.grid(True, which="both", alpha=0.2)

    # ---- B: the mechanism -- overlap falls exponentially with level ----
    n = np.arange(len(lv))
    for rs, c in [(0.15, "C2"), (0.25, "C0"), (0.4, "C3")]:
        O = np.sqrt(gen.overlap_masses(lv, source_size=rs))
        axB.semilogy(n, O / O[0], "o-", color=c, lw=1.5, ms=5,
                     label=f"source size {rs}  (×{gen.geometric_factor(O**2):.1f}/level)")
    axB.set_xlabel("internal level index $n$  (arithmetic ladder)")
    axB.set_ylabel(r"overlap $O_n/O_0$")
    axB.set_title("Overlap falls exponentially with level\n"
                  r"$\Rightarrow$ mass$\,\sim O_n^2$ is geometric", fontsize=11)
    axB.legend(fontsize=8.5, loc="upper right")
    axB.grid(True, which="both", alpha=0.2)

    fig.suptitle("Generations: an arithmetic ladder of internal levels becomes a "
                 "geometric mass hierarchy", fontsize=12.5, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "generations.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "generations.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'generations.pdf')}")


if __name__ == "__main__":
    main()
