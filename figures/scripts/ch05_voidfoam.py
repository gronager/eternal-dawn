#!/usr/bin/env python3
r"""The void as a nucleation-and-growth foam, and the average OGU size.

Left: a Johnson-Mehl tessellation -- OGUs nucleate in the void and grow at the
causal speed until they impinge, tiling it into a foam (impinged cells are curved
polygons; an isolated young OGU would be round). Right: the average OGU mass set by
the birth rate, M_avg ~ c^3 t_fill/2G with t_fill ~ (c/beta)^{1/4} -- frequent
births give many small universes, rare births few large ones.
Renders figures/pdf/void_foam.pdf, writes sims/output/void_foam.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import void_foam as vf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)
M_OUR = 9.2e52


def main() -> None:
    labels, areas, _ = vf.johnson_mehl_2d(n_seeds=90, grid=260, seed=3)
    real = areas[areas > 1e-4]

    lines = [
        "The void as a nucleation-and-growth foam: the average OGU size",
        "=" * 60,
        "  OGUs nucleate at rate beta and grow at c until they impinge (KJMA):",
        "    t_fill = (3/(pi beta c^3))^{1/4},  L ~ c t_fill,  M_avg ~ c^3 t_fill/2G.",
        "",
        "  birth rate beta needed for a given average OGU mass:",
        "    M_avg [kg]      beta [/m^3/s]     t_fill [s]      cell L [m]",
    ]
    for M in (1e54, 1e59, 1e65):
        b = vf.birth_rate_for_mass(M)
        lines.append(f"    {M:.0e}      {b:.2e}      {vf.fill_time(b):.2e}    "
                     f"{vf.char_size(b):.2e}")
    lines += [
        "",
        f"  Monte-Carlo foam (2D Johnson-Mehl): {len(real)} mature cells; cell-size",
        f"  spread std/mean = {real.std()/real.mean():.2f} (a broad but peaked",
        f"  distribution -- older cells larger, like the OGU age spread).",
        "",
        "READING: an individual OGU never runs out (infinite void, steady horizon",
        "growth), but the POPULATION does: OGUs tile the void and impinge, and that",
        "collective 'running out of nothing' fixes a characteristic cell size. The",
        "average OGU mass is then the birth rate read through the growth law,",
        "M_avg ~ c^3 (c/beta)^{1/4}/2G: frequent births -> many small universes, rare",
        "births -> few large ones. So 'how big is the average OGU' = 'how often is an",
        "OGU born' (Q12) -- the same single unknown the OGU mass and supraverse age",
        "reduce to. Pin beta (or measure one OGU size) and the whole foam scale follows.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "void_foam.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: the foam tessellation (random colours per cell).
    rng = np.random.default_rng(7)
    n = labels.max() + 1
    perm = rng.permutation(n)
    shuffled = perm[labels]
    axL.imshow(shuffled, origin="lower", extent=(0, 1, 0, 1),
               cmap="nipy_spectral", interpolation="nearest")
    axL.set_xticks([])
    axL.set_yticks([])
    axL.set_title("The condensed void as foam (Johnson-Mehl)\n"
                  "OGUs nucleate, grow at $c$, impinge -- cells are polygons, not "
                  "circles", fontsize=10.5)

    # Panel R: average OGU mass vs birth rate.
    betas = np.logspace(-150, -95, 300)
    M = np.array([vf.avg_ogu_mass(b) for b in betas])
    axR.plot(betas, M, "C0", lw=2.2)
    axR.axhspan(1e54, 1e65, color="C2", alpha=0.13)
    axR.text(betas[10], 3e58, "our OGU range\n($10^{54}$--$10^{65}$ kg)",
             fontsize=8, color="C2")
    for M0, lbl in [(1e54, "$\\sim$only child"), (1e65, "$\\sim10^{11}$ siblings")]:
        b0 = vf.birth_rate_for_mass(M0)
        axR.plot([b0], [M0], "ko", ms=5)
        axR.annotate(lbl, xy=(b0, M0), xytext=(6, -2), textcoords="offset points",
                     fontsize=7.5)
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"OG birth rate  $\beta$  [m$^{-3}$ s$^{-1}$]")
    axR.set_ylabel(r"average OGU mass  $M_{\rm avg}\sim c^3 t_{\rm fill}/2G$  [kg]")
    axR.set_title(r"Average OGU size $=$ birth rate via growth"
                  "\n"
                  r"$M_{\rm avg}\sim (c/\beta)^{1/4}\,c^2/2G$ (rarer births $\to$ bigger)",
                  fontsize=10.5)
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "void_foam.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
