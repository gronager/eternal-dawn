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


LOGRHO_C = -29.0       # the cosmic / critical density universes settle at (g cm^-3)


def logR_bh(logM):
    return logM + BH_C


def logR_co(logM):
    return CO_C - logM


A_GAP, TAU_DS, TAU_HK = 26.6, 12.0, 8.0   # TARDIS gap at the seed; de Sitter vs Hawking decay scales


def _envelope(logM):
    """The floor of the wedge: the upper of the two exclusion lines (membrane high-M, Compton low-M),
    smoothly rounded at the Planck apex where they meet."""
    k = 2.0
    return np.logaddexp(k * logR_bh(logM), k * logR_co(logM)) / k


def logR_in(logM):
    """The TARDIS inside view as a FULL ARC between the two exclusions. HIGH mass: asymptotes to the
    membrane (BH) line -- de Sitter, Omega -> 1, inside = outside (the OGU). LOW mass: asymptotes to the
    Compton line -- the Hawking/quantum limit, where an evaporating hole's 'inside' has shrunk to a
    single quantum (a torsiton), inside = outside again, the other way. The TARDIS bulge between is
    'bigger on the inside'; it peaks at the void seed, the smallest full universe. Both asymptotes meet
    at the Planck apex. Schematic."""
    logM = np.asarray(logM, float)
    bump = np.where(logM >= M_SEED,
                    A_GAP * np.exp(-(logM - M_SEED) / TAU_DS),     # de Sitter side, decays upward
                    A_GAP * np.exp(-(M_SEED - logM) / TAU_HK))     # Hawking side, decays downward
    out = _envelope(logM) + bump
    return out if out.ndim else float(out)


def main():
    fig, ax = plt.subplots(figsize=(9.8, 11.2))
    xlo, xhi, ylo, yhi = -36, 42, -38, 66

    # ---- forbidden regions + boundaries + Planck apex ------------------------------------------
    yy = np.linspace(ylo, yhi, 400)
    ax.fill_betweenx(yy, xlo, yy + BH_C, color="#7a3b3b", alpha=0.15, lw=0)
    xx = np.linspace(xlo, xhi, 400)
    ax.fill_between(xx, ylo, CO_C - xx, color="#5a5a5a", alpha=0.17, lw=0)
    ax.plot(yy + BH_C, yy, color="0.1", lw=1.7)                                   # membrane / OUTSIDE line
    ax.plot(xx, CO_C - xx, color="0.15", lw=1.4, ls=(0, (5, 2)))                  # Compton line
    ax.plot([R_PL], [M_PL], "s", color="black", ms=7)
    ax.annotate("Planck apex = smallest possible\nblack hole; the de Sitter and Hawking\n"
                "asymptotes of the TARDIS arc meet here", (R_PL, M_PL), xytext=(6, -34),
                textcoords="offset points", fontsize=8.3, color="0.2")

    # cosmic epoch lines (constant density, slope 3) -- our genesis epochs; NO GUT in ED
    def epoch_line(logrho, label, lx):
        rr = np.linspace(xlo, xhi, 60)
        mm = 3.0 * rr + logrho + 0.62
        ax.plot(rr, mm, color="0.5", lw=0.8, ls=(0, (4, 3)), alpha=0.7, zorder=1)
        ly = 3.0 * lx + logrho + 0.62
        ax.text(lx, ly, label, fontsize=7.8, color="0.4", rotation=72, rotation_mode="anchor",
                ha="left", va="bottom")
    epoch_line(15.0, "torsitonisation (mass-on)", -7.5)
    epoch_line(2.0, "nuclear (BBN)", -2.0)
    epoch_line(-12.0, "atomic (recomb)", 3.0)

    # ---- energy scales on the right (E = mc^2): the Planck/bounce scale, electroweak, and the CMB ---
    GEV = 23.75                                                # log10(1 g . c^2 / GeV)
    for elab, lgev, col in [(r"$E_P$ (the bounce scale)", 19.09, "0.4"),
                            (r"$E_{\rm EW}$", 2.0, "0.4")]:
        m = lgev - GEV
        ax.axhline(m, color=col, lw=0.6, ls=":", alpha=0.5)
        ax.text(xhi - 0.5, m + 0.5, elab, fontsize=8, color=col, ha="right")
    m_cmb = -12.6 - GEV                                        # CMB photon energy ~2.7 K
    ax.axhline(m_cmb, color="C5", lw=0.9, ls=":", alpha=0.8)
    ax.text(xlo + 1.5, m_cmb + 0.7, r"CMB $\approx$ the parent's Hawking radiation",
            fontsize=9, color="C5", fontweight="bold")
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

    # ---- the full Weltformel spectrum on the Compton limit: ALL 12 fermions (the torsiton rungs),
    #      the composite W/Z/H "mesons" (torsiton-antitorsiton pairs), and the nucleons p, n -------
    def lc(m):                                                 # Compton radius for a given log-mass
        return CO_C - m
    groups = [
        ("quarks", "o", "C1", [("u", -26.41), ("d", -26.08), ("s", -24.78),
                               ("c", -23.64), ("b", -23.13), ("t", -21.51)]),
        ("charged leptons", "*", "C4", [("e", -27.04), (r"$\mu$", -24.72), (r"$\tau$", -23.50)]),
        ("neutrinos", "v", "C9", [(r"$\nu_1$", -34.5), (r"$\nu_2$", -34.0), (r"$\nu_3$", -33.5)]),
        (r"$W,Z,H$ (composite)", "D", "C2", [("W", -21.84), ("Z", -21.79), ("H", -21.65)]),
    ]
    for name, mk, col, items in groups:
        for lab, m in items:
            ax.plot([lc(m)], [m], mk, color=col, ms=(11 if mk == "*" else 7.5), zorder=5)
            ax.annotate(lab, (lc(m), m), xytext=(-2, 4), textcoords="offset points", fontsize=7.4,
                        color=col, ha="right", fontweight="bold")
    for lab, m, r in [("p", -23.78, -13.05), ("n", -23.78, -12.78)]:                # composite nucleons
        ax.plot([r], [m], "s", color="C0", ms=6, zorder=5)
        ax.annotate(lab, (r, m), xytext=(3, -3), textcoords="offset points", fontsize=7.4, color="C0")

    # ---- the TARDIS: OUTSIDE (membrane line) vs INSIDE (cosmic, asymptoting to the line) ---------
    # the inside line ASYMPTOTES to the membrane: the TARDIS gap is enormous for a young void seed and
    # narrows up the line, the inside converging on the outside at the OGU (de Sitter, Omega -> 1).
    mline = np.linspace(M_SEED, M_OGU, 60)
    ax.plot([logR_in(m) for m in mline], mline, color="C0", lw=2.2, alpha=0.9, zorder=4)  # INSIDE line
    # the HAWKING tail: the same arc continued DOWN below the seed -- as a hole evaporates its inside
    # shrinks toward a single quantum, the curve peeling off the membrane to asymptote onto the Compton
    # (quantum) line at the fermion scale. The two asymptotes (de Sitter, Hawking) meet at the apex.
    mline_hk = np.linspace(-26.0, M_SEED, 60)
    ax.plot([logR_in(m) for m in mline_hk], mline_hk, color="C0", lw=1.5, ls=(0, (4, 2)),
            alpha=0.75, zorder=4)
    ax.annotate("the Hawking tail: evaporating, the inside\n"
                "shrinks to a single quantum (a torsiton) and\n"
                "the arc asymptotes to the Compton line --\n"
                "inside $=$ outside again, the quantum way",
                (logR_in(1.0), 1.0), xytext=(7, -2), textcoords="offset points",
                fontsize=7.8, color="C0")
    for m, col, big in [(M_SEED, "C2", False), (34.0, "0.3", False),
                        (56.0, "C3", False), (M_OGU, "C1", True)]:
        ro, ri = logR_bh(m), logR_in(m)
        ax.plot([ro, ri], [m, m], color=col, lw=1.0, ls=":", zorder=4)            # TARDIS connector
        ax.plot([ro], [m], "o", color=col, ms=(11 if big else 8), zorder=6)       # OUTSIDE (dense)
        ax.plot([ri], [m], "o", color=col, ms=(9 if big else 7), mfc="white", zorder=6)  # INSIDE (thin)

    # the smallest MASS-formed black hole: where torsitonisation (nuclear density) meets the BH line
    mb_x, mb_y = 6.1, 33.9
    ax.plot([mb_x], [mb_y], "P", color="0.1", ms=9, zorder=6)
    ax.annotate("smallest MASS-formed black hole\n($\\sim$ few $M_\\odot$, nuclear density) --\n"
                "below here a BH is compressed ENERGY\n(the tally), not compressed mass",
                (mb_x, mb_y), xytext=(10, -2), textcoords="offset points", fontsize=8.0, color="0.15")

    # the HAWKING WATERSHED: T_H(M) = T_CMB (~Moon mass, 4.5e25 g). Below it an UNFED hole net-
    # evaporates faster than it grows -- the fork in the journey. The void seed sits BELOW it: what
    # carries the seed UP across the membrane is the parent's feeding, not its own mass.
    m_eq = 25.65
    ax.plot([logR_bh(m_eq)], [m_eq], "X", color="#d2691e", ms=11, zorder=7)
    ax.annotate("Hawking watershed: $T_H\\!=\\!T_{\\rm CMB}$ ($\\sim$Moon mass)\n"
                "below here an UNFED hole net-evaporates (Hawking)\n"
                "instead of growing -- the seed is carried up by the\n"
                "parent's feeding, not by its own mass", (logR_bh(m_eq), m_eq),
                xytext=(10, 3), textcoords="offset points", fontsize=7.8, color="#d2691e")

    # the HAWKING RADIATION line: the wavelength of the EMITTED quanta (lambda ~ 16 R_s, Wien peak) --
    # PARALLEL to the membrane, ~1.2 dex above it. As a sub-watershed hole evaporates DOWN the membrane,
    # its emitted radiation tracks this line. (Distinct from the Compton/quantum tail of the inside arc.)
    HK_OFF = 1.2
    mh = np.linspace(M_PL, 34.0, 60)
    ax.plot(logR_bh(mh) + HK_OFF, mh, color="#d2691e", lw=1.3, ls=(0, (1, 2)), alpha=0.85, zorder=3)
    ax.text(logR_bh(8.0) + HK_OFF, 8.0, r"Hawking radiation ($\lambda\approx16\,R_s$)", fontsize=7.8,
            color="#d2691e", rotation=45, rotation_mode="anchor", ha="left", va="bottom")

    ax.annotate("void seed (mountain mass) --\ncompressed energy in the void", (logR_bh(M_SEED), M_SEED),
                xytext=(8, -16), textcoords="offset points", fontsize=8.4, color="C2", fontweight="bold")
    ax.annotate("a smaller black hole:\ntiny outside, a universe inside", (logR_bh(34), 34),
                xytext=(8, -22), textcoords="offset points", fontsize=8.0, color="0.3")
    ax.annotate("our universe -- the CRITICAL point\n(inside $=$ outside, $R_{\\rm Hubble}\\!=\\!R_s$)",
                (logR_bh(56), 56), xytext=(7, -20), textcoords="offset points", fontsize=8.6,
                color="C3", fontweight="bold")
    ax.annotate("OGU -- firstborn, largest; inside $\\to$ outside\n"
                "(de Sitter, $\\Omega\\!\\to\\!1$: the inside asymptotes to the line)", (logR_bh(M_OGU), M_OGU),
                xytext=(9, -6), textcoords="offset points", fontsize=9.0, color="C1", fontweight="bold")
    # the "tally" wedge: between the membrane and torsitonisation lines, below their meeting -- no mass
    ax.text(-1.0, 23.0, "the tally:\nno mass yet", fontsize=9, color="#6a3d9a", style="italic",
            ha="center", va="center", fontweight="bold")

    # ---- the inner-density JOURNEY: a field carrying each universe from cradle to grave ----------
    # born at rho_C (the no-mass tally density, on the critical line) -> expands, thinning -> dies at
    # near-void density (far right, under-critical). Not one line: a flow toward lower density.
    for m in (34.0, 41.0, 48.0, 55.0):
        x0 = max(logR_bh(m), 8.0) + 1.5
        ax.annotate("", xy=(x0 + 9.0, m - 3.0), xytext=(x0, m),
                    arrowprops=dict(arrowstyle="-|>", color="0.5", lw=1.2, alpha=0.5))
    ax.text(29.5, 25.0, "the inner-density journey:\nborn at $\\rho_C$ (cradle, no mass) $\\to$\n"
            "expand, thinning $\\to$ grave (near-void)", fontsize=8.2, color="0.4", va="top")

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

    # ---- the 12 fermions as a 4x3 table (towers x generations) in the open lower-right ----------
    tx = [30.0, 35.0, 40.0]                                   # generation columns I, II, III
    ty = [-5.0, -10.0, -15.0, -20.0]                          # tower rows: up, down, lepton, neutrino
    ax.text(17.0, 4.0, "the 12 fermions $=$ the torsiton generations", fontsize=10,
            color="0.12", fontweight="bold")
    ax.text(17.0, 1.0, "4 towers (charge/colour grip) $\\times$ 3 rungs (the bag-sharpness ladder)",
            fontsize=8.2, color="0.35")
    for j, g in enumerate(["gen I", "gen II", "gen III"]):
        ax.text(tx[j], -1.2, g, fontsize=8.6, color="0.3", ha="center", fontweight="bold")
    rows = [("up-type quark", "o", "C1", ["u", "c", "t"]),
            ("down-type quark", "o", "C1", ["d", "s", "b"]),
            ("charged lepton", "*", "C4", ["e", r"$\mu$", r"$\tau$"]),
            ("neutrino", "v", "C9", [r"$\nu_1$", r"$\nu_2$", r"$\nu_3$"])]
    for i, (rl, mk, col, names) in enumerate(rows):
        ax.text(17.0, ty[i], rl, fontsize=8.4, color=col, va="center")
        for j, nm in enumerate(names):
            ax.plot([tx[j] - 1.6], [ty[i]], mk, color=col, ms=(11 if mk == "*" else 7.5), zorder=5)
            ax.text(tx[j] - 0.6, ty[i], nm, fontsize=9, color=col, va="center", fontweight="bold")
    ax.text(17.0, -24.5, "their 15 masses (with the composite $W,Z,H$) are the torsiton result,\n"
            "created in the genesis cascade -- the same plane as the cosmology",
            fontsize=8.0, color="0.3", va="top")
    ax.plot([18], [-29.5], "o", color="0.3", ms=7)
    ax.text(19, -29.5, "outside (compact, dense)", fontsize=7.8, color="0.3", va="center")
    ax.plot([18], [-32.0], "o", color="0.3", ms=6, mfc="white")
    ax.text(19, -32.0, "inside (cosmic, thin) -- the TARDIS", fontsize=7.8, color="0.3", va="center")

    # ---- the TARDIS, made rigorous: a police-box-sized BLACK HOLE (R_s ~ 1 m) on the gravity line,
    #      with a baby universe on the inside line -- the only honest way to be bigger on the inside.
    m_t = 29.83                                                # R_s ~ 1 m (logR = 2) -> ~100 Earth masses
    ro_t, ri_t = logR_bh(m_t), logR_in(m_t)
    ax.plot([ro_t, ri_t], [m_t, m_t], color="#1f3a93", lw=1.0, ls=":", zorder=5)
    ax.plot([ro_t], [m_t], "s", color="#1f3a93", ms=8, zorder=6)
    ax.plot([ri_t], [m_t], "s", color="#1f3a93", ms=8, mfc="white", zorder=6)
    ax.annotate("TARDIS: a police-box\nblack hole ($R_s\\sim$1 m)", (ro_t, m_t), xytext=(-4, 9),
                textcoords="offset points", fontsize=7.8, color="#1f3a93", ha="right")
    ax.annotate("...a baby universe inside\n(the only honest way to be bigger on the inside)",
                (ri_t, m_t), xytext=(6, -2), textcoords="offset points", fontsize=7.8, color="#1f3a93")

    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, 73)
    ax.set_xlabel(r"$\log_{10}\,(\,$physical radius $/\,\mathrm{cm})$")
    ax.set_ylabel(r"$\log_{10}\,(\,$mass $/\,\mathrm{g})$")
    secax = ax.secondary_yaxis('right', functions=(lambda g: g + 23.75, lambda e: e - 23.75))
    secax.set_ylabel(r"$\log_{10}\,(\,$energy $/\,\mathrm{GeV})$    ($E=mc^2$)")
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
