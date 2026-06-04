#!/usr/bin/env python3
r"""Forces from soliton overlaps: colour channels (rigorous), and the confinement gap.

A: one-gluon-exchange colour-channel potentials V = <T1.T2> * alpha/r. The colour
   factors are RIGOROUS (forced by the framework's 3-valued label = SU(3) fundamental,
   identical to QCD); only the antisymmetric channels attract (q-qbar singlet -> meson,
   q-q antitriplet -> confined diquark); octet/sextet repel. ALL of them SCREEN
   (V -> 0). Dashed: the linear confinement V ~ sigma r that real quarks obey -- which
   the overlap/exchange picture does NOT produce, shown for contrast.
B: the residual force between two colour singlets from the overlap of their condensate
   fields -- an attractive Yukawa (the nuclear-force analogue), also SCREENED (-> 0).

The honest point: every overlap/exchange force here asymptotes to ZERO; CONFINEMENT
(~r) is non-perturbative (a flux tube), owed and undecided -- the deepest open question
of the strong sector. Renders figures/pdf/color_force.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import color_force as cf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    lines = [
        "Forces from overlaps: colour channels (rigorous) + the confinement gap",
        "=" * 70,
        f"  SU(3) validation: <T.T> fundamental = {cf.casimir_fundamental():.4f} (=4/3)",
        "  colour factors <T1.T2> (forced by the label; identical to QCD):",
    ]
    for key in ("qqbar_singlet", "qqbar_octet", "qq_antitriplet", "qq_sextet"):
        v = cf.color_factor(key)
        lines.append(f"    {cf.CHANNELS[key][3]:28s} {v:+.3f}  "
                     f"{'ATTRACT' if v < 0 else 'repel'}")
    lines += [
        "",
        "  free colour singlets only: q-qbar (2-body meson) or qqq (3-body baryon);",
        "  no free 2-quark state (diquark is antitriplet, not singlet).",
        "",
        "  THE CONFINEMENT GAP: every overlap/exchange force here SCREENS (V -> 0 at",
        "  large r). Real quarks are CONFINED (V ~ sigma r, never zero) -- a",
        "  non-perturbative flux tube the overlap picture does NOT produce. Whether the",
        "  gravity-torsion connection confines (~r) or only screens (->0) is the deepest",
        "  open question; if it only screens, quarks are not confined and the picture",
        "  fails here. Lattice-scale, undecided, honestly flagged.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "color_force.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    r = np.linspace(0.2, 6.0, 500)
    styles = {"qqbar_singlet": ("C3", "-", "q-q̄ singlet → meson (attract)"),
              "qq_antitriplet": ("C0", "-", "q-q antitriplet (confined diquark)"),
              "qqbar_octet": ("0.5", "--", "q-q̄ octet (repel)"),
              "qq_sextet": ("0.7", ":", "q-q sextet (repel)")}
    for key, (c, ls, lab) in styles.items():
        axA.plot(r, cf.oge_potential(r, key), color=c, ls=ls, lw=2.0, label=lab)
    axA.plot(r, -cf.linear_confinement(r, 0.25), "C1-.", lw=1.6,
             label="confinement $\\sim\\sigma r$ (real quarks; NOT from overlap)")
    axA.axhline(0, color="k", lw=0.8)
    axA.set_ylim(-2.0, 1.2)
    axA.set_xlabel(r"separation $r$ (soliton units)")
    axA.set_ylabel(r"$V(r)$")
    axA.set_title("Colour channels: factors RIGOROUS, but all SCREEN ($\\to0$)\n"
                  "confinement ($\\sim r$, dashed) is NOT produced -- owed", fontsize=10.5)
    axA.legend(fontsize=8, loc="lower left")
    axA.grid(True, alpha=0.2)

    axB.plot(r, cf.residual_yukawa(r, g=2.5, m_sigma=1.0), "C2", lw=2.2,
             label="residual (overlap σ-exchange)")
    axB.axhline(0, color="k", lw=0.8)
    axB.set_xlabel(r"separation between two singlets $r$")
    axB.set_ylabel(r"residual $V(r)$")
    axB.set_title("Residual singlet-singlet force (nuclear analogue)\n"
                  "attractive, SCREENED ($\\to0$) -- correct for the nuclear force",
                  fontsize=10.5)
    axB.legend(fontsize=8.5, loc="lower right")
    axB.grid(True, alpha=0.2)

    fig.suptitle("Forces as overlap channels: the colour structure is forced; but "
                 "overlaps SCREEN -- confinement ($\\sim r$) is the owed gap", fontsize=11.5,
                 y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "color_force.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "color_force.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'color_force.pdf')}")


if __name__ == "__main__":
    main()
