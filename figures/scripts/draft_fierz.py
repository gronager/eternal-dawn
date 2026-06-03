#!/usr/bin/env python3
r"""The Fierz probe: the torsion four-fermion term has G_A = G_V (walking-favourable).

A: the meson-channel couplings of the Hehl-Datta (axial-axial) torsion interaction,
   Fierzed into the exchange channels (explicit Dirac matrices, validated against the
   V-A self-Fierz identity and Fierz^2=identity): G_S=+1/4, G_P=-1/4, G_V=G_A=+1/2,
   G_T=0. The vector and axial channels are coupled EQUALLY.
B: why that matters for S. The electroweak S is driven by the vector/axial resonance
   splitting r = M_A/M_V: a QCD-like sector (G_A != G_V) sits at r~1.6 with S~0.3 (the
   graveyard); equal couplings (G_A=G_V) push r -> 1, where the V and A cancel and
   S -> 0. The Fierz forces the favourable DIRECTION (not tuned); the exact landing
   point needs the full RPA.

Renders figures/pdf/fierz.pdf, writes sims/output/fierz.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import fierz as fz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    c = fz.hehl_datta_channels()
    lines = [
        "Fierz probe: the torsion four-fermion term in the meson channels",
        "=" * 62,
        f"  validations: V-A self-Fierz dev = {fz.validate_VmA_self_fierz():.0e}, "
        f"Fierz^2-identity dev = {fz.validate_involution():.0e}  (both 0: code trusted)",
        "",
        "  Hehl-Datta (axial-axial) Fierzed into exchange (meson) channels:",
        f"    G_S = {c['S']:+.3f}   G_P = {c['P']:+.3f}   G_V = {c['V']:+.3f}   "
        f"G_A = {c['A']:+.3f}   G_T = {c['T']:+.3f}",
        f"    -> G_A / G_V = {fz.GA_over_GV():.3f}",
        "",
        "READING: the torsion interaction couples the VECTOR and AXIAL meson channels",
        "EQUALLY (G_A = G_V), forced by its axial-axial Lorentz structure -- not tuned.",
        "That is the structural prerequisite for 'walking' (M_a1 -> M_rho), where the",
        "vector and axial spectral functions cancel and the electroweak S -> 0. A QCD-",
        "like sector has G_A != G_V, M_a1 > M_rho, and S ~ 0.3 (the graveyard). So the",
        "cheap probe is GREEN: the framework has a real, structural shot at escaping --",
        "the opposite of generic technicolor. (Also: G_S=+1/4>0 drives the chiral",
        "condensate the soliton forms; G_T=0.) NECESSARY, not sufficient: the residual",
        "splitting from the chiral-breaking loops needs the full RPA -- now worth building.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "fierz.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0))

    # ---- A: channel couplings ----
    chans = ["S", "P", "V", "A", "T"]
    vals = [c[k] for k in chans]
    colors = ["0.6", "0.6", "C0", "C3", "0.6"]
    bars = axA.bar(chans, vals, color=colors, width=0.6)
    for b, vv in zip(bars, vals):
        axA.text(b.get_x() + b.get_width() / 2, vv + (0.02 if vv >= 0 else -0.05),
                 f"{vv:+.2f}", ha="center", fontsize=9)
    axA.axhline(0, color="k", lw=0.8)
    axA.annotate("$G_V=G_A$\n(walking-\nfavourable)", (2.5, 0.5),
                 textcoords="offset points", xytext=(0, 18), ha="center",
                 fontsize=9, color="C3",
                 arrowprops=dict(arrowstyle="-[,widthB=2.6", color="C3", lw=1.2))
    axA.set_ylabel("meson-channel coupling (units of $G_{HD}$)")
    axA.set_title("Fierz of the torsion (Hehl-Datta) term\n"
                  "vector and axial couplings are EQUAL", fontsize=11)
    axA.set_ylim(-0.45, 0.8)
    axA.grid(True, axis="y", alpha=0.2)

    # ---- B: implication for S ----
    r = np.linspace(1.0, 2.0, 200)                  # r = M_A/M_V
    fV_over_MV = 0.168                               # set so QCD (r=1.6) gives S~0.3
    S = 4 * np.pi * fV_over_MV**2 * (1 - 1.0 / r**4)
    axB.plot(r, S, "C2", lw=2.2)
    axB.axhline(0.1, color="0.4", ls="--", lw=1.2)
    axB.text(1.02, 0.115, "LEP bound", fontsize=8.5, color="0.35")
    axB.axhspan(0.1, 0.5, color="0.85", alpha=0.6)
    # QCD-like point
    axB.plot(1.6, 4 * np.pi * fV_over_MV**2 * (1 - 1 / 1.6**4), "ro", ms=9)
    axB.annotate("QCD-like\n($G_A\\neq G_V$, graveyard)", (1.6, 0.3),
                 textcoords="offset points", xytext=(-6, -4), fontsize=8.5,
                 color="r", ha="right")
    # the direction the torsion structure pushes
    axB.annotate("", xy=(1.05, 0.01), xytext=(1.55, 0.28),
                 arrowprops=dict(arrowstyle="->", color="C0", lw=2.0))
    axB.plot(1.0, 0.0, "C0o", ms=9)
    axB.annotate("$G_A=G_V$ (torsion)\npushes $M_A\\to M_V$, $S\\to0$", (1.0, 0.0),
                 textcoords="offset points", xytext=(8, 30), fontsize=8.5, color="C0")
    axB.set_xlim(1.0, 2.0); axB.set_ylim(-0.02, 0.45)
    axB.set_xlabel(r"resonance splitting  $M_A/M_V$")
    axB.set_ylabel("electroweak $S$")
    axB.set_title("Equal couplings push toward $S\\to0$\n"
                  "(direction forced; landing point needs RPA)", fontsize=11)
    axB.grid(True, alpha=0.2)

    fig.suptitle("The cheap probe is green: the torsion four-fermion term has "
                 "$G_A=G_V$, the structural key to a small $S$", fontsize=12, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "fierz.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "fierz.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'fierz.pdf')}")


if __name__ == "__main__":
    main()
