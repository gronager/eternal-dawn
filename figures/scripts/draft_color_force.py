#!/usr/bin/env python3
r"""Forces from soliton overlaps: the colour channels and the residual nuclear force.

A: the colour-channel potentials V = <T1.T2> x Cornell(r). The colour factors are
   RIGOROUS (forced by the framework's 3-valued antisymmetric label = SU(3)
   fundamental); the radial shape is modelled (the string tension is lattice-scale).
   Only the antisymmetric channels attract: q-qbar singlet (mesons) and q-q
   antitriplet (the confined diquark); octet and sextet repel. So free states are
   colour singlets -- q-qbar (2-body meson) or qqq (3-body baryon) -- and there is no
   free 2-quark state, exactly as observed.
B: the residual force between two colour SINGLETS from the overlap of their condensate
   fields -- an attractive Yukawa, the nuclear-force analogue. This is the 'force from
   overlaps' the soliton picture delivers directly.

Renders figures/pdf/color_force.pdf, writes sims/output/color_force.txt.
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
        "Forces from soliton overlaps: colour channels (rigorous) + residual (overlap)",
        "=" * 76,
        f"  SU(3) validation: <T.T> fundamental = {cf.casimir_fundamental():.4f} (= 4/3)",
        "",
        "  colour factors <T1.T2> (forced by the 3-valued label; identical to QCD):",
    ]
    for key in ("qqbar_singlet", "qqbar_octet", "qq_antitriplet", "qq_sextet"):
        _, _, _, lab = cf.CHANNELS[key]
        v = cf.color_factor(key)
        lines.append(f"    {lab:28s} {v:+.3f}   {'ATTRACT' if v < 0 else 'repel'}")
    lines += [
        "",
        "  => free colour singlets only: q-qbar (2-body MESON) or qqq (3-body BARYON);",
        "     no free 2-quark state (diquark is antitriplet, not singlet).",
        "",
        "STATUS: the colour-channel structure is RIGOROUS and forced by the framework's",
        "Pauli-required label. The radial shape (Cornell) is MODELLED -- the string",
        "tension is lattice-scale, not computed. The residual singlet-singlet force is",
        "the overlap-derived nuclear-force analogue (sigma-exchange Yukawa).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "color_force.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.2))

    # ---- A: colour-channel potentials ----
    r = np.linspace(0.12, 2.2, 400)
    styles = {"qqbar_singlet": ("C3", "-", "q-q̄ singlet → meson"),
              "qq_antitriplet": ("C0", "-", "q-q antitriplet (confined diquark)"),
              "qqbar_octet": ("0.5", "--", "q-q̄ octet (repulsive)"),
              "qq_sextet": ("0.7", ":", "q-q sextet (repulsive)")}
    for key, (c, ls, lab) in styles.items():
        axA.plot(r, cf.channel_potential(r, key), color=c, ls=ls, lw=2.0, label=lab)
    axA.axhline(0, color="k", lw=0.8)
    axA.set_ylim(-2.2, 2.2)
    axA.set_xlabel(r"separation $r$ (soliton units)")
    axA.set_ylabel(r"$V(r)=\langle T_1\!\cdot\! T_2\rangle\times$Cornell$(r)$")
    axA.set_title("Colour-channel forces\n(factors RIGOROUS; shape modelled)", fontsize=11)
    axA.legend(fontsize=8.5, loc="upper left")
    axA.grid(True, alpha=0.2)
    axA.text(1.2, -1.7, "only antisymmetric\nchannels bind", fontsize=8.5, color="0.4")

    # ---- B: residual force from overlap (nuclear-force analogue) ----
    rr = np.linspace(0.25, 4.0, 400)
    axB.plot(rr, cf.residual_yukawa(rr, g=2.0, m_sigma=1.0), "C2", lw=2.2)
    axB.axhline(0, color="k", lw=0.8)
    axB.set_xlabel(r"separation between two singlets $r$")
    axB.set_ylabel(r"residual $V(r)$ (overlap)")
    axB.set_title("Residual force from the overlap of two singlets\n"
                  "(the nuclear-force analogue, $\\sigma$-exchange)", fontsize=11)
    axB.grid(True, alpha=0.2)
    axB.annotate("attractive,\nshort-range", (1.2, cf.residual_yukawa(1.2, 2.0, 1.0)),
                 textcoords="offset points", xytext=(20, -6), fontsize=8.5, color="C2")

    fig.suptitle("Forces as overlap channels: the colour structure is forced by the "
                 "label; the residual is the nuclear force", fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "color_force.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "color_force.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'color_force.pdf')}")


if __name__ == "__main__":
    main()
