#!/usr/bin/env python3
r"""The supraverse census: viable universes by type, and our footprint.

The wallpaper math tabulates viable observers along three independent axes --
generation depth (geometric), chirality (50/50 by CPT), and family (clean vs
baryon-rich) -- and locates where we sit. Left: the generation tower with our cell
marked. Right: the (chirality x family) breakdown at BHU1, showing our footprint holds
~40% of all viable observers (a typical, Copernican place). Renders
figures/pdf/census.pdf, writes sims/output/census.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import census as cs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

EPS = 0.1
F_CLEAN = 0.9


def main():
    fp = cs.our_footprint(EPS, F_CLEAN)
    n, pgen = cs.generation_pmf(EPS, n_max=6)

    lines = [
        "The supraverse census: viable universes by type, and our footprint",
        "=" * 66,
        "  three independent axes (each from the wallpaper math):",
        "   generation n : geometric P(n|n>=1)=(1-eps)eps^{n-1}  (anchor + branching)",
        "   chirality     : 50/50 matter/antimatter (CPT) -- the halves are identical",
        "   family        : clean (eta~ppb, fair-sample progenitor) vs baryon-rich",
        "",
        "  generation tower (eps=0.1):",
    ]
    for ni, pn in zip(n, pgen):
        bar = "#" * int(round(pn * 50))
        lines.append(f"    BHU{int(ni)}:  {pn:6.4f}  {bar}")
    lines += [
        "",
        f"  OUR FOOTPRINT: BHU{fp.generation}, {fp.chirality}, {fp.family} family",
        f"    = {fp.fraction:.3f} of all viable observers (~{fp.fraction*100:.0f}%)",
        f"    typical (dominant cell)? {cs.footprint_is_typical(EPS, F_CLEAN)}",
        "",
        "  observable tags of our cell:",
        "    generation : BHU1 (DM+DE anchor n>=1, ~90% are n=1)",
        "    chirality  : matter (our half; antimatter twin identical by CPT)",
        "    family     : clean (eta~6e-10 -> fair-sample/horizon-scale progenitor)",
        "    size       : ~a Hubble mass (Nariai-capped ~ M_vis)",
        "    spin       : low/quiet typical -- but LOCAL, tags the progenitor hole",
        "                 NOT the generation (a null axis does NOT make us BHU2)",
        "",
        "READING: our derived tags (BHU1, matter, clean, ~Hubble mass) land in the",
        "DOMINANT census cell -- ~40% of all viable observers share them. That is the",
        "Copernican footprint check: we are where most observers are, not a rare",
        "outlier. The one tag that is NOT a generation marker is spin: it is the local",
        "progenitor hole's, re-drawn each generation, so a shared-axis null leaves our",
        "BHU1 placement untouched -- it only says our progenitor hole was quiet.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "census.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: generation tower.
    ns = [int(x) for x in n[:5]]
    ps = pgen[:5]
    colors = ["C0" if k == 1 else "0.6" for k in ns]
    axL.bar(ns, ps, color=colors, alpha=0.85, width=0.65)
    for k, p in zip(ns, ps):
        axL.text(k, p + 0.015, f"{p:.3f}", ha="center", fontsize=8.5)
    axL.text(1, 0.9 * 0.5, "we are here\n(BHU1)", ha="center", fontsize=9,
             color="white", fontweight="bold")
    axL.set_xlabel(r"generation depth  $n$  (BHU$_n$)")
    axL.set_ylabel(r"$P(\mathrm{BHU}_n \mid n\geq1)$")
    axL.set_title("Generation census: geometric, BHU1-dominated\n"
                  r"$P(n)=(1-\epsilon)\epsilon^{n-1}$, $\epsilon=0.1$", fontsize=11)
    axL.set_xticks(ns)
    axL.set_ylim(0, 1.0)
    axL.grid(True, axis="y", alpha=0.2)

    # Panel R: BHU1 breakdown by chirality x family (a 2x2 mosaic).
    axR.set_xlim(0, 1)
    axR.set_ylim(0, 1)
    axR.set_aspect("equal")
    axR.axis("off")
    axR.set_title("Within BHU1: chirality $\\times$ family\n"
                  "our cell (blue) holds $\\sim$40% of viable observers", fontsize=11)
    p1 = pgen[0]                                   # BHU1 share
    # columns = chirality (matter|antimatter), rows = family (clean top, baryon-rich bot)
    cells = [
        (0.0, 0.0, 0.5, 0.9, "C0", "matter\nclean\n(US)", True),
        (0.5, 0.0, 0.5, 0.9, "0.55", "antimatter\nclean\n(CPT twin)", False),
        (0.0, 0.9, 0.5, 0.1, "0.78", "matter\nbaryon-rich", False),
        (0.5, 0.9, 0.5, 0.1, "0.82", "antimatter\nbaryon-rich", False),
    ]
    from matplotlib.patches import Rectangle
    for x, y, w, h, c, lab, ours in cells:
        axR.add_patch(Rectangle((x, y), w, h, facecolor=c, edgecolor="white",
                                lw=2))
        axR.text(x + w / 2, y + h / 2, lab, ha="center", va="center",
                 fontsize=8.5 if h > 0.2 else 7,
                 color="white" if (ours or c < "0.6") else "0.1",
                 fontweight="bold" if ours else "normal")
    axR.text(0.5, -0.06, "(BHU1 = 90% of all observers; this square is that 90%)",
             ha="center", fontsize=8, color="0.4")

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "census.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "census.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'census.pdf')}")


if __name__ == "__main__":
    main()
