#!/usr/bin/env python3
r"""The lumpiness ladder: large-scale inhomogeneity accumulating across membranes.

Left: lumpiness vs generation -- void (0) -> OGU (minimal, the LambdaCDM smooth-start
case) -> BHU1 -> BHU2 -> ... converging to a fixed point. The biggest jump is
OGU->BHU1, where the parent's lumpiness first appears. Right: observed ultra-large
structures all exceed the homogeneity scale, sitting above the smooth-start (OGU)
floor -- evidence we are a BHU (n>=1), the excess being the inherited rung.
Renders figures/pdf/lumpiness.pdf, writes sims/output/lumpiness.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import lumpiness as lp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    L = lp.lumpiness_ladder()
    Lstar = lp.fixed_point()
    jump = lp.biggest_jump_index(L)

    lines = [
        "The lumpiness ladder: inhomogeneity across Cartasis membranes",
        "=" * 62,
        "  L_n^2 = a^2 + (b L_{n-1})^2  (intrinsic + transferred parent lumpiness)",
        "",
        "  generation     lumpiness L",
        f"    void          {0.0:.2f}   (homogeneous, unstable to fluctuations)",
        f"    OGU (n=0)     {L[0]:.2f}   (no parent -> minimal; = LambdaCDM smooth start)",
    ]
    for n in range(1, len(L)):
        lines.append(f"    BHU{n}          {L[n]:.2f}")
    lines += [
        f"    fixed point   {Lstar:.2f}   (deep generations saturate, b<1)",
        "",
        f"  biggest jump: {'OGU->BHU1' if jump==0 else f'BHU{jump}->BHU{jump+1}'}"
        f" (+{np.diff(L)[jump]:.2f}) -- the parent's lumpiness first appears.",
        "",
        "  observed ultra-large structures vs the homogeneity scale (~260 Mpc):",
    ]
    for name, mpc in sorted(lp.ULTRA_LARGE.items(), key=lambda x: -x[1]):
        flag = "ABOVE" if lp.exceeds_homogeneity(mpc) else "below"
        lines.append(f"    {name:20s} {mpc:6.0f} Mpc   {flag} homogeneity scale")
    lines += [
        "",
        "READING: LambdaCDM assumes a smooth/singular start and builds structure only",
        "internally -- in SCT terms it predicts an OGU, the minimal-lumpiness rung. But",
        "we observe MORE: ultra-large structures (Giant Arc ~1 Gly, Big Ring ~0.4 Gpc,",
        "Sloan Great Wall, KBC void) exceed the homogeneity scale and are hard to",
        "assemble internally in the time available. In SCT that excess is the INHERITED",
        "rung of the ladder -- the parent's lumpiness projected through the membrane --",
        "so it is evidence we are a BHU (n>=1), not an OGU. Lumpiness is thus yet",
        "another generation-depth indicator, alongside dark matter and dark energy.",
        "The model is illustrative (a, b not yet derived from membrane projection); the",
        "robust, falsifiable content is the DIRECTION: monotonic excess large-scale",
        "power for n>=1, growing across membranes toward a fixed point.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "lumpiness.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: the ladder.
    xs = list(range(len(L)))
    labels = ["void", "OGU"] + [f"BHU{n}" for n in range(1, len(L) - 1)]
    yvals = [0.0] + list(L)
    xx = list(range(len(yvals)))
    colors = ["0.7", "0.5"] + ["C0"] * (len(yvals) - 2)
    axL.plot(xx, yvals, "-", color="0.6", lw=1.4, zorder=1)
    axL.scatter(xx, yvals, c=colors, s=80, zorder=3)
    axL.axhline(Lstar, color="C3", ls="--", lw=1.2)
    axL.text(len(yvals) - 1.2, Lstar + 0.04, f"fixed point $L^*={Lstar:.1f}$",
             fontsize=8.5, color="C3", ha="right")
    # mark the big OGU->BHU1 jump and "us"
    axL.annotate("parent's lumpiness\nfirst appears\n(biggest jump)",
                 xy=(2, yvals[2]), xytext=(2.6, 0.55), fontsize=8,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.scatter([2], [yvals[2]], s=160, facecolor="none", edgecolor="C0", lw=2,
                zorder=4)
    axL.text(2, yvals[2] + 0.12, "us (BHU1)", fontsize=8.5, color="C0", ha="center")
    axL.axvspan(-0.4, 0.4, color="0.85", alpha=0.5)
    axL.text(0, 0.06, "homogeneous\nvoid", fontsize=8, color="0.4", ha="center")
    axL.text(1, yvals[1] - 0.16, "OGU\n($\\Lambda$CDM\nsmooth start)", fontsize=7.5,
             color="0.35", ha="center")
    axL.set_xticks(xx)
    axL.set_xticklabels(["void", "OGU"] + [f"BHU{n}" for n in range(1, len(L))],
                        rotation=30, fontsize=8)
    axL.set_ylabel("large-scale lumpiness  $L$  (excess power)")
    axL.set_title("Lumpiness climbs a ladder across membranes\n"
                  "void $\\to$ OGU $\\to$ BHU1 $\\to$ ... $\\to L^*$", fontsize=11)
    axL.set_ylim(0, Lstar * 1.25)
    axL.grid(True, axis="y", alpha=0.2)

    # Panel R: observed ultra-large structures vs homogeneity scale.
    items = sorted(lp.ULTRA_LARGE.items(), key=lambda x: x[1])
    names = [k for k, _ in items]
    vals = [v for _, v in items]
    cols = ["C0" if lp.exceeds_homogeneity(v) else "0.6" for v in vals]
    axR.barh(names, vals, color=cols, alpha=0.85)
    axR.axvline(lp.HOMOGENEITY_MPC, color="C3", ls="--", lw=1.5)
    axR.text(lp.HOMOGENEITY_MPC + 8, 0.1, "homogeneity\nscale ~260 Mpc\n"
             "(smooth-start / OGU floor)", fontsize=8, color="C3")
    axR.set_xlabel("comoving extent  [Mpc]")
    axR.set_title("Observed ultra-large structure exceeds the floor\n"
                  "the excess = inherited lumpiness $\\Rightarrow$ we are a BHU",
                  fontsize=11)
    axR.grid(True, axis="x", alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "lumpiness.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "lumpiness.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'lumpiness.pdf')}")


if __name__ == "__main__":
    main()
