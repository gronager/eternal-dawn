#!/usr/bin/env python3
r"""Chapter 9 (Seeing almost Nothing): the Eternal Dawn mass--radius diagram.

The cosmic mass--radius plane reworked for ED, built around the TARDIS. Every object lives on the
diagonal of existence between two exclusions: the TORSION--GRAVITY MEMBRANE (the black-hole line
R = 2GM/c^2; cross it and torsion bounces you into a child universe) and the COMPTON limit (quantum
uncertainty), meeting at the Planck apex.

THE TARDIS. A black hole is TWO points at one mass: the compact OUTSIDE (on the membrane line, high
density rho ~ 1/M^2) and the cosmic INSIDE (a whole universe, far larger, far thinner). The gap is
enormous for a mountain-mass void seed -- microscopic outside, a universe inside -- and SHRINKS up
the line, the inside and outside converging at the OGU (the firstborn, largest before it evaporates),
which is so vast it is its own interior. Beyond it: the scaleless void.

The genesis cascade (mass-on -> confinement -> baryons -> nuclei) happens INSIDE the young universe
and CREATES the particles on the diagonal of existence -- so it is marked there, at the particle
scale, not on the universe line. Renders figures/pdf/ed_massradius.pdf. Schematic, cgs.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

BH_C = -27.83          # log R_bh = log M + BH_C     (R = 2GM/c^2)
CO_C = -37.46          # log R_co = -log M + CO_C     (R = hbar/Mc)
M_PL, R_PL = -4.74, -32.79
M_SEED, M_OGU = 16.0, 63.0


def logR_bh(logM):
    return logM + BH_C


def logR_in(logM):
    """The cosmic INSIDE view: outside R_s plus a TARDIS gap that is huge at the seed and closes
    (converges onto the membrane line) at the OGU. Schematic."""
    gap = 40.0 * (M_OGU - logM) / (M_OGU - M_SEED)        # 40 dex at the seed -> 0 at the OGU
    return logR_bh(logM) + gap


def main():
    fig, ax = plt.subplots(figsize=(9.8, 10.0))
    xlo, xhi, ylo, yhi = -36, 42, -30, 66

    # ---- forbidden regions + boundaries + Planck apex ------------------------------------------
    yy = np.linspace(ylo, yhi, 400)
    ax.fill_betweenx(yy, xlo, yy + BH_C, color="#7a3b3b", alpha=0.15, lw=0)
    xx = np.linspace(xlo, xhi, 400)
    ax.fill_between(xx, ylo, CO_C - xx, color="#5a5a5a", alpha=0.17, lw=0)
    ax.plot(yy + BH_C, yy, color="0.1", lw=1.7)                                   # membrane / OUTSIDE line
    ax.plot(xx, CO_C - xx, color="0.15", lw=1.4, ls=(0, (5, 2)))                  # Compton line
    ax.plot([R_PL], [M_PL], "s", color="black", ms=7)
    ax.annotate("Planck apex =\nsmallest possible\nblack hole", (R_PL, M_PL), xytext=(6, -30),
                textcoords="offset points", fontsize=8.3, color="0.2")

    # cosmic epoch lines (constant density, slope 3) -- our genesis epochs; NO GUT in ED
    def epoch_line(logrho, label, lx):
        rr = np.linspace(xlo, xhi, 60)
        mm = 3.0 * rr + logrho + 0.62
        ax.plot(rr, mm, color="0.5", lw=0.8, ls=(0, (4, 3)), alpha=0.7, zorder=1)
        ly = 3.0 * lx + logrho + 0.62
        ax.text(lx, ly, label, fontsize=7.8, color="0.4", rotation=72, rotation_mode="anchor",
                ha="left", va="bottom")
    epoch_line(15.0, "QGP", -7.5)
    epoch_line(2.0, "nuclear (BBN)", -2.0)
    epoch_line(-12.0, "atomic (recomb)", 3.0)
    ax.text(-27, 49, "torsion--gravity\nmembrane\n(cross it $\\to$ a child\nuniverse bounces out)",
            fontsize=9.5, color="#7a3b3b", ha="center", va="center", style="italic")
    ax.text(-18, -23, "quantum uncertainty", rotation=-43, fontsize=10, color="0.3",
            ha="center", va="center", style="italic")

    # ---- diagonal of existence: representative objects (the INSIDE of our universe) ------------
    objs = [(-8, -23, "atom"), (-5, -15, "virus"), (2, 5, "human"), (8.8, 27.8, "Earth"),
            (10.8, 33.3, "Sun"), (9, 33.3, "WD"), (6, 33.6, "NS"), (22.7, 45.3, "Milky Way"),
            (25, 48, "clusters")]
    ax.plot([o[0] for o in objs], [o[1] for o in objs], "o", color="C0", ms=4, alpha=0.8, zorder=3)
    for x, y, lab in objs:
        ax.annotate(lab, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7.5, color="C0")

    # ---- the charged-lepton spectrum on the Compton limit: OUR mass spectrum (the torsiton rungs)
    leptons = [(-10.42, -27.04, "e"), (-12.71, -24.73, r"$\mu$"), (-13.96, -23.50, r"$\tau$")]
    ax.plot([l[0] for l in leptons], [l[1] for l in leptons], "*", color="C4", ms=11, zorder=5)
    for x, y, lab in leptons:
        ax.annotate(lab, (x, y), xytext=(-3, 4), textcoords="offset points", fontsize=9,
                    color="C4", ha="right", fontweight="bold")
    ax.annotate("the charged-lepton spectrum $e,\\mu,\\tau$\n"
                "(the torsiton generations -- our mass result)\ncreated in the genesis cascade",
                (-12.7, -24.7), xytext=(2.0, -16.0), textcoords="data", fontsize=8.2, color="C4",
                ha="left", va="top", arrowprops=dict(arrowstyle="->", color="C4", lw=0.9))

    # ---- the TARDIS: OUTSIDE (membrane line) vs INSIDE (cosmic), across the range ---------------
    ax.plot([logR_in(m) for m in np.linspace(M_SEED, M_OGU, 40)],
            np.linspace(M_SEED, M_OGU, 40), color="C0", lw=2.2, alpha=0.9, zorder=4)   # INSIDE line
    pairs = [(M_SEED, "void seed (mountain mass)\nmicroscopic outside, a universe inside", "C2"),
             (34.0, "a smaller black hole", "0.3"),
             (56.0, "our universe", "C3"),
             (M_OGU, "OGU (firstborn, largest\nbefore it evaporates)", "C1")]
    for m, lab, col in pairs:
        ro, ri = logR_bh(m), logR_in(m)
        ax.plot([ro, ri], [m, m], color=col, lw=1.0, ls=":", zorder=4)            # TARDIS connector
        ax.plot([ro], [m], "o", color=col, ms=8, zorder=6)                        # OUTSIDE (dense)
        if m < M_OGU - 0.1:
            ax.plot([ri], [m], "o", color=col, ms=7, mfc="white", zorder=6)       # INSIDE (thin)
    # the seed (a diamond) + our universe labels
    ax.annotate("void seed (mountain mass)", (logR_bh(M_SEED), M_SEED), xytext=(8, -14),
                textcoords="offset points", fontsize=8.6, color="C2", fontweight="bold")
    ax.annotate("a smaller black hole:\ntiny outside, a universe inside", (logR_bh(34), 34),
                xytext=(8, -20), textcoords="offset points", fontsize=8.0, color="0.3")
    ax.annotate("our universe -- OUTSIDE\n(a black hole in the parent)", (logR_bh(56), 56),
                xytext=(7, -20), textcoords="offset points", fontsize=8.6, color="C3", fontweight="bold")
    ax.annotate("OGU -- inside $\\approx$ outside\n(it is its own interior)", (logR_bh(M_OGU), M_OGU),
                xytext=(9, -6), textcoords="offset points", fontsize=9.2, color="C1", fontweight="bold")

    # ---- the DARK SECTOR: child universes are still fed, so the OUTSIDE mass CLIMBS the line ----
    for m in (M_SEED, 34.0, 56.0):
        ro = logR_bh(m)
        ax.add_patch(FancyArrowPatch((ro - 0.5, m + 0.3), (ro + 1.2, m + 2.0), arrowstyle="-|>",
                                     mutation_scale=12, color="#c000c0", lw=2.0, zorder=7))
    ax.annotate("still fed by the parent: the OUTSIDE mass climbs the line --\n"
                "dark matter (the mass arriving) $+$ dark energy (its inflow rate) $=$ the ongoing dawn\n"
                "(the OGU is parentless: no inflow, no dark sector -- it only evaporates)",
                (logR_bh(56), 56), xytext=(-2.0, 62.0), textcoords="data", fontsize=8.0,
                color="#c000c0", ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#c000c0", lw=0.9))

    # the two density labels along the two lines
    ax.text(logR_bh(61) - 5.5, 61 + 1.2, r"OUTSIDE: $\rho \propto 1/M^2$ (compact, dense)",
            fontsize=8.6, color="0.2", rotation=45, rotation_mode="anchor")
    ax.text(logR_in(40) - 1.0, 40 + 0.5, "INSIDE: cosmic, thin\n(bigger on the inside)",
            fontsize=8.6, color="C0", rotation=45, rotation_mode="anchor", ha="left")

    # the void
    ax.text(2.0, 70.0, "the void", fontsize=14, color="0.2", style="italic", fontweight="bold")
    ax.text(2.0, 67.4, "infinite, scaleless -- the nested tower floats in it", fontsize=8.6, color="0.4")
    ax.add_patch(FancyArrowPatch((logR_bh(M_OGU) + 0.4, M_OGU + 0.8), (14.0, 69.0), arrowstyle="->",
                                 mutation_scale=13, color="0.5", lw=1.1))

    leg = [Line2D([], [], color="0.1", lw=2, label="OUTSIDE: universes as black holes (the membrane)"),
           Line2D([], [], color="C0", lw=2, label="INSIDE: the cosmos each contains (the TARDIS view)"),
           Line2D([], [], color="0.3", marker="o", lw=0, label="outside point (compact, dense)"),
           Line2D([], [], color="0.3", marker="o", mfc="white", lw=0, label="inside point (cosmic, thin)")]
    ax.legend(handles=leg, loc="lower right", fontsize=8.2, framealpha=0.92)

    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, 73)
    ax.set_xlabel(r"$\log_{10}\,(\,$physical radius $/\,\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}\,(\,$mass $/\,\mathrm{g})$")
    ax.set_title("The Eternal Dawn mass--radius diagram\n"
                 "the TARDIS: outside (compact) and inside (a universe) converge at the OGU; beyond, the void",
                 fontsize=11.5)
    ax.grid(alpha=0.12)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "ed_massradius.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")
    for m in (M_SEED, 34, 56, M_OGU):
        print(f"  logM={m:.0f}: outside logR={logR_bh(m):.1f}, inside logR={logR_in(m):.1f}, "
              f"gap={logR_in(m)-logR_bh(m):.0f} dex")


if __name__ == "__main__":
    main()
