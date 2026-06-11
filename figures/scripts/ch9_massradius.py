#!/usr/bin/env python3
r"""Chapter 9 (Seeing almost Nothing): the Eternal Dawn mass--radius diagram.

The cosmic mass--radius plane reworked for ED. Every object lives on the diagonal of existence between
two exclusions: the TORSION--GRAVITY MEMBRANE (the black-hole line R = 2GM/c^2; cross it and torsion
bounces you into a child universe) and the COMPTON limit (quantum uncertainty), meeting at the Planck
apex.

THE ED STORY IS ONE ASCENDING LINE. A universe expands while staying near-critical (horizon mass ~
Schwarzschild mass), so its whole life is a climb UP the membrane line:
  * the void seed (a mountain-mass black hole -- the cheapest way in) sits low on the line;
  * particle production at the bounce grows the mass up the line through the genesis cascade
    (mass-on -> confinement -> baryons -> nuclei) -- our own fermion/baryon creation lines;
  * the domination eras (radiation -> matter -> dark energy) are SEGMENTS of that climb, with matter =
    the parent's mass through the membrane and dark energy = the parent's ongoing inflow;
  * "now" is our Hubble point; the line continues through the parent to the OGU (firstborn, largest
    before it evaporates), and beyond it lies the scaleless void.
TARDIS: each black hole is TWO points at one mass -- the compact OUTSIDE (on the line) and the cosmic
INSIDE (a whole universe, to the right) -- which converge onto the line at the cosmic scale.
Renders figures/pdf/ed_massradius.pdf. Schematic, orders of magnitude, cgs.
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

BH_C = -27.83          # log R_bh = log M + BH_C    (R = 2GM/c^2)
CO_C = -37.46          # log R_co = -log M + CO_C    (R = hbar/Mc)
M_PL, R_PL = -4.74, -32.79


def logR_bh(logM):
    return logM + BH_C


def main():
    fig, ax = plt.subplots(figsize=(9.6, 10.2))
    xlo, xhi, ylo, yhi = -36, 42, -30, 66

    # ---- forbidden regions ---------------------------------------------------------------------
    yy = np.linspace(ylo, yhi, 400)
    ax.fill_betweenx(yy, xlo, yy + BH_C, color="#7a3b3b", alpha=0.16, lw=0)        # membrane side
    xx = np.linspace(xlo, xhi, 400)
    ax.fill_between(xx, ylo, CO_C - xx, color="#5a5a5a", alpha=0.18, lw=0)         # quantum side

    # ---- domination eras as SEGMENTS of the climb up the membrane line (labelled in the legend) -
    eras = [(16.0, 45.0, "#d98a8a", r"$\Omega_r$  radiation"),
            (45.0, 52.0, "#8aa0d9", r"$\Omega_m$  matter (parent's mass through the membrane)"),
            (52.0, 63.0, "#9a9a9a", r"$\Omega_\Lambda$  dark energy (parent's ongoing inflow)")]
    for m0, m1, col, lab in eras:
        mm = np.linspace(m0, m1, 30)
        ax.plot(logR_bh(mm), mm, color=col, lw=12, alpha=0.55, solid_capstyle="butt", zorder=2)

    # ---- the two boundary lines + Planck apex --------------------------------------------------
    ax.plot(yy + BH_C, yy, color="0.1", lw=1.7)                                    # membrane / BH line
    ax.plot(xx, CO_C - xx, color="0.15", lw=1.5, ls=(0, (5, 2)))                   # Compton line
    ax.plot([R_PL], [M_PL], "o", color="black", ms=6)
    ax.annotate("Planck apex", (R_PL, M_PL), xytext=(6, -14), textcoords="offset points",
                fontsize=8.5, color="0.2")
    ax.text(-26, 50, "torsion--gravity\nmembrane\n(cross it $\\to$ a child\nuniverse bounces out)",
            fontsize=9.5, color="#7a3b3b", ha="center", va="center", style="italic")
    ax.text(-19, -23, "quantum uncertainty", rotation=-43, fontsize=10, color="0.3",
            ha="center", va="center", style="italic")

    # ---- the genesis cascade (our fermion/baryon creation lines), low on the climb -------------
    stages = [(20.0, r"$\rho_C$ membrane"), (27.0, "mass-on"), (33.0, "confinement"),
              (39.0, "baryons$+$leptons"), (44.5, "nuclei (BBN)")]
    for m, lab in stages:
        R = logR_bh(m)
        ax.plot([R - 1.4, R + 1.4], [m + 1.4, m - 1.4], color="0.25", lw=0.8)      # short crossing tick
        ax.annotate(lab, (R - 1.4, m + 1.4), xytext=(-3, 1), textcoords="offset points",
                    fontsize=7.6, color="0.25", ha="right")

    # ---- diagonal of existence: representative objects -----------------------------------------
    objs = [(-12.6, -27, "e"), (-13, -23.8, "p"), (-8, -23, "atom"), (-5, -15, "virus"),
            (2, 5, "human"), (8.8, 27.8, "Earth"), (10.8, 33.3, "Sun"), (9, 33.3, "WD"),
            (6, 33.6, "NS"), (22.7, 45.3, "Milky Way"), (25, 48, "clusters")]
    ax.plot([o[0] for o in objs], [o[1] for o in objs], "o", color="C0", ms=4, alpha=0.8, zorder=3)
    for x, y, lab in objs:
        ax.annotate(lab, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7.5, color="C0")

    # ---- the void seed (mountain mass): the cheapest way in ------------------------------------
    M_seed, R_seed = 16.0, logR_bh(16.0)
    ax.plot([R_seed], [M_seed], "D", color="C2", ms=9, zorder=6)
    ax.annotate("void seed\n(mountain mass --\nthe cheapest way in)", (R_seed, M_seed),
                xytext=(10, -6), textcoords="offset points", fontsize=8.5, color="C2", fontweight="bold")
    # particle production grows the mass UP the line
    ax.annotate("particle production\nat the bounce grows\nthe mass up the line",
                (logR_bh(30), 30), xytext=(14, -8), textcoords="offset points", fontsize=8,
                color="C2", arrowprops=dict(arrowstyle="->", color="C2", lw=0.9))

    # ---- our universe (TARDIS two dots), parent, OGU, void -------------------------------------
    M_us, R_us = 56.0, logR_bh(56.0)
    M_par, R_par = 59.0, logR_bh(59.0)
    M_ogu, R_ogu = 63.0, logR_bh(63.0)
    R_in = R_us + 6.0                                   # the "inside" (cosmic) view: bigger on the inside
    ax.plot([R_us, R_in], [M_us, M_us], color="C3", lw=1.0, ls=":", zorder=5)
    ax.plot([R_us], [M_us], "o", color="C3", ms=9, zorder=6)
    ax.plot([R_in], [M_us], "o", color="C3", ms=7, mfc="white", zorder=6)
    ax.annotate("our universe: OUTSIDE\n(a black hole in the parent)", (R_us, M_us), xytext=(6, -22),
                textcoords="offset points", fontsize=8.6, color="C3", fontweight="bold")
    ax.annotate("INSIDE (our cosmos --\nbigger on the inside)", (R_in, M_us), xytext=(6, 2),
                textcoords="offset points", fontsize=8.2, color="C3")
    ax.plot([R_par], [M_par], "o", color="C3", ms=8, zorder=6)
    ax.annotate("parent", (R_par, M_par), xytext=(10, -2), textcoords="offset points",
                fontsize=8.5, color="C3")
    ax.plot([R_ogu], [M_ogu], "o", color="C1", ms=11, zorder=6)
    ax.annotate("OGU -- firstborn,\nlargest before it dies", (R_ogu, M_ogu), xytext=(9, -10),
                textcoords="offset points", fontsize=9.5, color="C1", fontweight="bold")

    ax.text(logR_bh(61.0) - 5.5, 61.0 + 1.2, r"$\rho \propto 1/M^2$", fontsize=9, color="0.2",
            rotation=45, rotation_mode="anchor")

    # the void
    ax.text(0.0, 70.0, "the void", fontsize=14, color="0.2", style="italic", fontweight="bold")
    ax.text(0.0, 67.4, "infinite, scaleless -- the nested tower floats in it", fontsize=8.6, color="0.4")
    ax.add_patch(FancyArrowPatch((R_ogu + 0.4, M_ogu + 0.8), (12.0, 69.0), arrowstyle="->",
                                 mutation_scale=13, color="0.5", lw=1.1))

    # legend: the domination eras (segments of the climb), TARDIS, and the seed
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=c, lw=9, alpha=0.6, label=lab) for (_, _, c, lab) in eras]
    handles += [Line2D([], [], color="C3", marker="o", lw=0, label="black hole: OUTSIDE (on the line)"),
                Line2D([], [], color="C3", marker="o", mfc="white", lw=0,
                       label="same one: INSIDE (a universe) -- TARDIS, converging at the cosmic scale"),
                Line2D([], [], color="C2", marker="D", lw=0, label="void seed (mountain mass)")]
    ax.legend(handles=handles, loc="lower right", fontsize=8.0, framealpha=0.92,
              title="the one ascending line", title_fontsize=8.5)

    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, 73)
    ax.set_xlabel(r"$\log_{10}\,(\,$physical radius $/\,\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}\,(\,$mass $/\,\mathrm{g})$")
    ax.set_title("The Eternal Dawn mass--radius diagram\n"
                 "one ascending line: void seed $\\to$ genesis $\\to$ us $\\to$ OGU, then the void",
                 fontsize=12)
    ax.grid(alpha=0.12)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "ed_massradius.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
