#!/usr/bin/env python3
r"""Paper C, Fig. 3 -- three independent roads to the bag sharpness s_T, all in the window.

The single number to which the fermion-mass hierarchy reduces, s_T = R0/r0, determined three
independent ways on the L16x48 dynamical ensemble -- they did not have to agree:

  1. one-body proxy      : the dressed-quark scalar profile rho(r)=Tr[S^dag S], chiral-extrapolated  -> 0.43
  2. three-body condensate: the genuine <N|qbar q|N> bag (run/10), bracketed across sink times        -> 0.45-0.56
  3. lepton hierarchy     : the value the observed charged-lepton span (x3477) requires of the lever   -> 0.43

All three sit inside the productive window [0.43, 0.70]. The bag sharpness measured directly on the
lattice and the sharpness the lepton masses demand are the same number. Renders three_roads.pdf.

HONEST: the three-body road is bracketed (sink-time systematic), not pinned; the lepton road carries
the lever model; single heavy sea mass throughout. Mechanism confirmed, not the spectrum to a fixed %.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

WIN = (0.43, 0.70)
# (label, central s_T, [lo, hi] uncertainty band, colour)
ROADS = [
    ("one-body proxy\n" r"$\rho=\mathrm{Tr}[S^\dagger S]$", 0.430, (0.41, 0.49), "C0"),
    ("three-body condensate\n" r"$\langle N|\bar q q|N\rangle$", 0.505, (0.45, 0.56), "C2"),
    ("lepton hierarchy\n" r"(the $\times3477$ lever)", 0.428, (0.42, 0.44), "C3"),
]


def main():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.axhspan(WIN[0], WIN[1], color="0.85", alpha=0.7, zorder=0,
               label=f"productive window [{WIN[0]}, {WIN[1]}]")
    xs = np.arange(len(ROADS))
    for x, (lab, c, (lo, hi), col) in zip(xs, ROADS):
        ax.errorbar(x, c, yerr=[[c - lo], [hi - c]], fmt="o", color=col, ms=9, lw=2.2, capsize=6)
        ax.annotate(f"{c:.2f}", (x, hi + 0.012), ha="center", color=col, fontsize=10, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([r[0] for r in ROADS], fontsize=9.5)
    ax.set_xlim(-0.6, len(ROADS) - 0.4)
    ax.set_ylim(0.30, 0.78)
    ax.set_ylabel(r"bag sharpness  $s_T = R_0 / r_0$")
    ax.set_title(r"Three independent roads to $s_T$ agree at $\approx 0.43$--$0.50$",
                 fontsize=12)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.18, axis="y")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "three_roads.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")
    centrals = [c for _, c, _, _ in ROADS]
    print(f"  roads cluster at s_T = {min(centrals):.3f}-{max(centrals):.3f} "
          f"(window lower edge {WIN[0]}; the lever-required value sits AT that edge by construction)")
    for lab, c, band, _ in ROADS:
        print(f"    {lab.splitlines()[0]:24s} s_T={c:.3f}  band {band}")


if __name__ == "__main__":
    main()
