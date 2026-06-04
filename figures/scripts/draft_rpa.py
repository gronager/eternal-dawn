#!/usr/bin/env python3
r"""RPA-dressed S: the torsion structure's advantage grows in the walking limit.

A: S vs Lambda/M (the walking knob -- larger = less chiral breaking) for the torsion
   couplings (G_A=G_V, from the Fierz) and a QCD-like sector (G_A=0.5 G_V). As the
   sector walks, the equal couplings keep S BOUNDED while the QCD-like S blows up.
B: the ratio S_torsion/S_QCD falls from ~1 (strong breaking) toward ~0.3 (deep
   walking): the equal-coupling advantage is real and grows.

HONEST: this RPA is one-sided (vector-resonance enhancement, which RAISES S; not the
full axial/Weinberg-sum-rule catch-up, which LOWERS it), so the ABSOLUTE values are
overestimates and never drop below the leading ~0.16. It establishes the RELATIVE
advantage (torsion << QCD-like), not the absolute escape S<0.1 -- which needs the full
chiral RPA. Renders figures/pdf/rpa.pdf, writes sims/output/rpa.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import rpa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    Lams = np.array([3, 5, 10, 30, 100, 300, 1000.0])
    res = [rpa.compare(L) for L in Lams]
    St = np.array([r["S_torsion"] for r in res])
    Sq = np.array([r["S_qcd"] for r in res])
    ratio = np.array([r["ratio"] for r in res])

    lines = [
        "RPA-dressed S: torsion (G_A=G_V) vs QCD-like (G_A=0.5 G_V)",
        "=" * 58,
        f"{'Lam/M':>7} {'PiA/PiV':>9} {'S_torsion':>11} {'S_qcd':>9} {'ratio':>7}",
    ]
    for r in res:
        lines.append(f"{r['Lam_over_M']:>7.0f} {r['PiA_over_PiV']:>9.3f} "
                     f"{r['S_torsion']:>11.3f} {r['S_qcd']:>9.3f} {r['ratio']:>7.3f}")
    lines += [
        "",
        "READING: in the walking limit (large Lambda/M) the loops equalise and the",
        "torsion's EQUAL couplings keep S bounded while the QCD-like S blows up -- the",
        "ratio falls to ~0.3 (torsion ~3x smaller S). The Fierz direction is confirmed",
        "by an independent calculation: the torsion structure is on the right side of S.",
        "",
        "HONEST LIMIT: this RPA is ONE-SIDED -- it has the vector enhancement (raises S)",
        "but not the full axial/Weinberg catch-up (lowers S), so the ABSOLUTE values are",
        "overestimates and never fall below the leading ~0.16. It shows the RELATIVE",
        "advantage, not the absolute escape S<0.1. That still needs the full chiral RPA",
        "(likely lattice). S remains the owed make-or-break; the torsion is favourably",
        "placed, not proven to win.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "rpa.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0))

    axA.loglog(Lams, Sq, "ro-", lw=1.8, ms=6, label=r"QCD-like ($G_A=0.5\,G_V$)")
    axA.loglog(Lams, St, "C0o-", lw=1.8, ms=6, label=r"torsion ($G_A=G_V$, Fierz)")
    axA.axhline(3 / (6 * np.pi), color="0.5", ls=":", lw=1.0)
    axA.text(3, 3 / (6 * np.pi) * 1.1, "leading $N_c/6\\pi$", fontsize=8, color="0.4")
    axA.set_xlabel(r"$\Lambda/M$  (walking knob: larger = less chiral breaking)")
    axA.set_ylabel("RPA-dressed $S$ (overestimate; relative use)")
    axA.set_title("Equal couplings keep $S$ bounded;\nQCD-like blows up as it walks",
                  fontsize=11)
    axA.legend(fontsize=9, loc="upper left")
    axA.grid(True, which="both", alpha=0.2)

    axB.semilogx(Lams, ratio, "C2o-", lw=2.0, ms=6)
    axB.axhline(1.0, color="0.5", ls="--", lw=1.0)
    axB.text(3, 1.02, "no advantage", fontsize=8, color="0.4")
    axB.set_ylim(0, 1.1)
    axB.set_xlabel(r"$\Lambda/M$ (walking)")
    axB.set_ylabel(r"$S_{\rm torsion}/S_{\rm QCD\text{-}like}$")
    axB.set_title("The equal-coupling advantage\ngrows as the sector walks", fontsize=11)
    axB.grid(True, which="both", alpha=0.2)

    fig.suptitle("RPA: the torsion's $G_A=G_V$ gives a real, growing advantage on $S$ "
                 "(relative; absolute escape still owed)", fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "rpa.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "rpa.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'rpa.pdf')}")


if __name__ == "__main__":
    main()
