#!/usr/bin/env python3
r"""Chapter 9 (Seeing almost Nothing): the Eternal Dawn mass--radius diagram.

The cosmic mass--radius plane -- every object lives on the narrow "diagonal of existence" between
two exclusions: the BLACK-HOLE line (R = 2GM/c^2, "forbidden by gravity" above it) and the COMPTON
limit (R = hbar/Mc, "quantum uncertainty" below it), which pinch shut at the Planck apex.

THE ED TWIST: the Hubble point is not the edge. Our whole universe sits ON the black-hole line -- a
black hole from outside -- so the line CONTINUES up-right, through the parent, to the OGU (the
firstborn, the largest universe before it evaporates). All universes are points on this one line,
identical at birth (the universal bounce density rho_C at the membrane) and differing only in how
far up they slid (rho ~ 1/M^2: bigger = thinner). Beyond the OGU is the void -- the infinite,
scaleless substrate the whole nested tower floats in. Renders figures/pdf/ed_massradius.pdf.

Schematic (orders of magnitude), in cgs: log10(R/cm) on x, log10(M/g) on y.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

# the two boundaries, in log10 cgs (R = 2GM/c^2 ; lambda_C = hbar/Mc), meeting at the Planck point
BH_C = -27.83          # log R_bh = log M + BH_C   (R = 2GM/c^2)
CO_C = -37.46          # log R_co = -log M + CO_C   (R = hbar/Mc)
M_PL, R_PL = -4.74, -32.79   # Planck apex (log m_P [g], log l_P [cm])


def logR_bh(logM):
    return logM + BH_C


def main():
    fig, ax = plt.subplots(figsize=(9.2, 9.8))

    xlo, xhi, ylo, yhi = -36, 42, -30, 66

    # ---- the forbidden regions ----------------------------------------------------------------
    yy = np.linspace(ylo, yhi, 400)
    x_bh = yy - (-BH_C)                                  # x on the BH line for mass yy: x = logM + BH_C
    x_bh = yy + BH_C
    # "forbidden by gravity" = inside the horizon = LEFT of the BH line (smaller R than R_s)
    ax.fill_betweenx(yy, xlo, x_bh, color="#7a3b3b", alpha=0.25, lw=0)
    # "quantum uncertainty" = below the Compton line
    xx = np.linspace(xlo, xhi, 400)
    y_co = -(xx) + CO_C                                  # Compton: logR = -logM + CO_C -> logM = CO_C - logR
    ax.fill_between(xx, ylo, y_co, color="#5a5a5a", alpha=0.22, lw=0)

    # ---- the two boundary lines + Planck apex -------------------------------------------------
    ax.plot(yy + BH_C, yy, color="0.15", lw=1.6)                       # BH line
    ax.plot(xx, CO_C - xx, color="0.15", lw=1.6, ls=(0, (5, 2)))       # Compton line
    ax.plot([R_PL], [M_PL], "o", color="black", ms=7)
    ax.annotate("Planck apex\n(quantum gravity)", (R_PL, M_PL), xytext=(6, -22),
                textcoords="offset points", fontsize=8.5, color="0.2")

    ax.text(-21, 52, "forbidden by gravity\n(inside the horizon)", rotation=43, fontsize=10.5,
            color="#7a3b3b", ha="center", va="center", style="italic")
    ax.text(-19, -23, "quantum uncertainty", rotation=-43, fontsize=10.5, color="0.3",
            ha="center", va="center", style="italic")

    # ---- the diagonal of existence: representative objects ------------------------------------
    objs = [
        (-12.6, -27.0, "e"), (-13.0, -23.8, "p"), (-8.0, -23.0, "atom"),
        (-5.0, -15.0, "virus"), (2.0, 5.0, "human"), (8.8, 27.8, "Earth"),
        (10.8, 33.3, "Sun"), (9.0, 33.3, "WD"), (6.0, 33.6, "NS"),
        (22.7, 45.3, "Milky Way"), (25.0, 48.0, "clusters"),
    ]
    ox = [o[0] for o in objs]; oy = [o[1] for o in objs]
    ax.plot(ox, oy, "o", color="C0", ms=4.5, alpha=0.85)
    for x, y, lab in objs:
        ax.annotate(lab, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=8, color="C0")

    # ---- ED: the universe line, our node, the OGU, the void ----------------------------------
    M_us, R_us = 56.0, logR_bh(56.0)                    # our universe, ON the BH line
    M_par, R_par = 59.0, logR_bh(59.0)                  # the parent
    M_ogu, R_ogu = 63.0, logR_bh(63.0)                  # the OGU (firstborn, largest before it dies)

    # the universe line continuing up-right from us to the OGU (heavier on the SAME BH line)
    ax.plot([R_us, R_ogu], [M_us, M_ogu], color="C3", lw=3.0, solid_capstyle="round", zorder=5)
    ax.plot([R_us, R_par], [M_us, M_par], "o", color="C3", ms=8, zorder=6)
    ax.plot([R_ogu], [M_ogu], "o", color="C1", ms=11, zorder=6)
    ax.annotate("our universe\n(a black hole from outside)", (R_us, M_us), xytext=(12, -20),
                textcoords="offset points", fontsize=9, color="C3", fontweight="bold")
    ax.annotate("parent", (R_par, M_par), xytext=(12, -3), textcoords="offset points",
                fontsize=8.5, color="C3")
    ax.annotate("OGU -- firstborn,\nlargest before it dies", (R_ogu, M_ogu), xytext=(9, -12),
                textcoords="offset points", fontsize=9.5, color="C1", fontweight="bold")

    # the density gradient along the line (placed along it, mid-upper, clear of the corner)
    ax.text(logR_bh(50.0) - 6.5, 50.0 + 1.5, r"$\rho \propto 1/M^2$ (bigger $=$ thinner)",
            fontsize=9.5, color="0.2", rotation=45, rotation_mode="anchor")

    # one consolidated note, in the open lower-right wedge (clear of the objects)
    ax.annotate("all universes are points on this one line:\n"
                r"identical at birth (one density $\rho_C$ at the membrane)," "\n"
                "sliding apart only in scale. Our Hubble point is a\n"
                "node on the parent's line -- and so on, up to the OGU.",
                (R_us, M_us), xytext=(16.0, 10.0), textcoords="data",
                fontsize=8.8, color="C3", ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="C3", lw=0.9))

    # the void: beyond the OGU, the scaleless substrate (top-centre, clear of the corner)
    ax.text(2.0, 70.2, "the void", fontsize=14, color="0.2", style="italic", fontweight="bold")
    ax.text(2.0, 67.6, "infinite, scaleless --\nthe nested tower floats in it",
            fontsize=8.8, color="0.4")
    arr = FancyArrowPatch((R_ogu + 0.4, M_ogu + 0.8), (12.0, 69.0), arrowstyle="->",
                          mutation_scale=13, color="0.5", lw=1.1)
    ax.add_patch(arr)

    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, 73)
    ax.set_xlabel(r"$\log_{10}\,(\,$physical radius $/\,\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}\,(\,$mass $/\,\mathrm{g})$")
    ax.set_title("The Eternal Dawn mass--radius diagram\n"
                 "the universe line runs through us to the OGU; beyond it, the void", fontsize=12)
    ax.grid(alpha=0.12)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "ed_massradius.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  our universe on the BH line: (logR={R_us:.1f}, logM={M_us:.0f}); "
          f"OGU at (logR={R_ogu:.1f}, logM={M_ogu:.0f}); rho ~ 1/M^2 along the line")


if __name__ == "__main__":
    main()
