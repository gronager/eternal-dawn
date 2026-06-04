#!/usr/bin/env python3
r"""The ab-initio spectrum from G, hbar, c alone (the Weltformel computation).

The 3x4 fermion matrix and W/Z/H computed with the three constants of Nature as the ONLY
dimensionful inputs -- no Yukawa, no v, no f_pi, no fitted well. Because the four-fermion
coupling IS gravity, the scale is the Planck mass, so every entry comes out ~10^17 too
heavy. Two honest gaps are shown:
  A: the 3x4 matrix (ab initio, GeV) -- all near Planck; the STRUCTURE (towers, generations,
     neutrino lightest) is there, the SCALE is Planck.
  B: ab initio vs observed (log-log) -- the horizontal offset is the hierarchy (the density
     question, L5); the compression is the generation spread the simple well underestimates.

Renders figures/pdf/ab_initio_spectrum.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import ab_initio_spectrum as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

TCOLOR = {"up-quark": "C3", "down-quark": "C2", "charged-lepton": "C0", "neutrino": "C4"}


def main():
    fm = sp.fermion_matrix_GeV()
    bm = sp.boson_masses_GeV()
    M_Pl = fm["M_Pl_GeV"]
    gap = sp.hierarchy_gap()

    lines = [
        "The ab-initio spectrum from G, hbar, c alone (Weltformel)",
        "=" * 60,
        f"  INPUT: G, hbar, c  ->  M_Pl c^2 = {M_Pl:.3e} GeV  (the ONLY scale)",
        f"  overlaps (no fit): {np.round(fm['overlaps'], 4)}  "
        f"gen-factor {np.round(fm['gen_factor'], 3)}",
        f"  handles (no fit): " + ", ".join(f"{k}={v:.4f}" for k, v in fm["handles"].items()),
        "",
        "  AB INITIO 3x4 FERMION MATRIX (GeV):",
        f"  {'tower':16s}{'gen I':>13s}{'gen II':>13s}{'gen III':>13s}",
    ]
    for t in sp.TOWERS:
        row = fm["matrix_GeV"][t]
        lines.append(f"  {t:16s}" + "".join(f"{v:13.3e}" for v in row)
                     + f"   ({','.join(sp.GEN_LABELS[t])})")
    lines += ["", "  OBSERVED (GeV):"]
    for t in sp.TOWERS:
        lines.append(f"  {t:16s}" + "".join(f"{v:13.3e}" for v in sp.OBSERVED_GEV[t]))
    lines += [
        "",
        f"  AB INITIO BOSONS (GeV): W={bm['W']:.3e}  Z={bm['Z']:.3e}  H={bm['H']:.3e}",
        f"    (observed W=80.4, Z=91.2, H=125.2)",
        f"  HIERARCHY GAP (the one number not produced): M_Pl/m_top = {gap:.2e}",
        "",
        "READING: this is the strict Weltformel -- one nonlinear Dirac (four-fermion) soliton,",
        "inputs G, hbar, c ONLY. The matrix is real and ab initio. Because the four-fermion",
        "coupling IS gravity, the scale is the PLANCK mass: every entry is ~10^17 too heavy.",
        "The STRUCTURE comes out (3 generations x 4 towers; quarks > leptons; neutrino the",
        "lightest, forced to zero in the singlet limit). Two honest gaps remain: the SCALE",
        "(Planck vs electroweak = the hierarchy = the density question L5) and the generation",
        "spread (~3 here vs ~10^4 observed -- the simple well; the steep suppression needs the",
        "real chiral/walking dynamics = lattice). Off by a lot, but turned out by the crank.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "ab_initio_spectrum.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.6))

    # ---- A: the 3x4 matrix as a grid of GeV values ----
    axA.set_xlim(0, 3); axA.set_ylim(0, 4)
    for j, t in enumerate(reversed(sp.TOWERS)):
        row = fm["matrix_GeV"][t]
        for i in range(3):
            axA.add_patch(plt.Rectangle((i, j), 1, 1, facecolor=TCOLOR[t], alpha=0.18,
                                        edgecolor="k", lw=0.6))
            axA.text(i + 0.5, j + 0.62, sp.GEN_LABELS[t][i], ha="center", fontsize=10,
                     fontweight="bold", color=TCOLOR[t])
            axA.text(i + 0.5, j + 0.32, f"{row[i]:.2e}", ha="center", fontsize=8, color="0.2")
    axA.set_xticks([0.5, 1.5, 2.5]); axA.set_xticklabels(["gen I", "gen II", "gen III"])
    axA.set_yticks([0.5, 1.5, 2.5, 3.5])
    axA.set_yticklabels([t.replace("-", "\n") for t in reversed(sp.TOWERS)], fontsize=8)
    axA.set_title(f"Ab-initio 3x4 fermion matrix (GeV)\nscale = $M_{{Pl}}$ = {M_Pl:.2e} GeV "
                  "(from $G,\\hbar,c$)", fontsize=11)
    axA.tick_params(length=0)

    # ---- B: ab initio vs observed (log-log) ----
    for t in sp.TOWERS:
        pred = fm["matrix_GeV"][t]
        obs = sp.OBSERVED_GEV[t]
        axB.scatter(obs, pred, s=60, color=TCOLOR[t], edgecolor="k", lw=0.5, zorder=3,
                    label=t)
    # bosons
    axB.scatter([80.4, 91.2, 125.2], [bm["W"], bm["Z"], bm["H"]], marker="s", s=70,
                color="C1", edgecolor="k", lw=0.5, zorder=3, label="W, Z, H")
    lo, hi = 1e-12, 1e20
    axB.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.5, label="if scale were right")
    axB.axhline(M_Pl, color="0.6", ls=":", lw=1.0)
    axB.text(1e-11, M_Pl * 1.5, "Planck scale (ab-initio output)", fontsize=8, color="0.4")
    axB.annotate("", xy=(1e2, 1e3), xytext=(1e2, M_Pl),
                 arrowprops=dict(arrowstyle="<->", color="C6", lw=1.4))
    axB.text(1.6e2, 1e10, f"hierarchy gap\n$\\sim${gap:.0e}\n(= density question, L5)",
             fontsize=8.5, color="C6")
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlim(1e-12, 1e4); axB.set_ylim(1e15, 1e20)
    axB.set_xlabel("observed mass (GeV)")
    axB.set_ylabel("ab-initio mass (GeV) -- Planck-scaled")
    axB.set_title("Ab initio vs observed\nstructure right, scale = Planck (off by the "
                  "hierarchy)", fontsize=11)
    axB.legend(fontsize=7.5, loc="lower left", ncol=2)
    axB.grid(True, which="both", alpha=0.15)

    fig.suptitle("The Weltformel computation: the whole spectrum from $G,\\hbar,c$ alone "
                 "-- structure ab initio, absolute scale = Planck (the hierarchy is owed)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(PDF_DIR, "ab_initio_spectrum.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "ab_initio_spectrum.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'ab_initio_spectrum.pdf')}")


if __name__ == "__main__":
    main()
