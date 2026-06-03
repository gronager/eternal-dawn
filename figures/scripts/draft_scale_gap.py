#!/usr/bin/env python3
r"""The scale gap as dimensional transmutation of the spin-gravity coupling.

A: Lambda/M_Pl = exp(-2pi/(b0 g_T^2)) -- the exponential collapse of the scale. The observed
   electroweak/Planck ratio (~10^-17) is hit by an ORDINARY coupling g_T^2 ~ O(0.01-0.1), not
   a tuning. The running coupling alpha_T(mu): weak at M_Pl, strong at Lambda (asympt. free).
B: the Weltformel 3x4 matrix BEFORE (Planck) and AFTER transmutation (Lambda) -- the heavy
   end drops from ~10^18 GeV into the observed 10-170 GeV ballpark. The scale gap closes;
   what remains is the generation spread (the substrate-overlap refinement, separate).

Renders figures/pdf/scale_gap.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import scale_gap as sg
from cartasis_sims import ab_initio_spectrum as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

TCOLOR = {"up-quark": "C3", "down-quark": "C2", "charged-lepton": "C0", "neutrino": "C4"}


def main():
    s = sg.scale_summary()
    M_Pl = s["M_Pl_GeV"]
    Lam = s["target_scale_GeV"]

    lines = [
        "The scale gap as dimensional transmutation of the spin-gravity coupling",
        "=" * 70,
        f"  M_Pl = {M_Pl:.3e} GeV   ->   target Lambda (top) = {Lam:.0f} GeV",
        f"  hierarchy = exp(-{s['ln_hierarchy']:.1f}) ~ 10^-17  (NOT a tuning -- exp of a coupling)",
        f"  coupling needed: g_T^2 = {s['g_T^2_needed']:.4f}, alpha_T = {s['alpha_T_needed']:.4f}"
        f"  (b0=7) -- ordinary",
        "",
        "  g_T^2 (and alpha_T) that hit the electroweak scale, for a range of b0:",
    ]
    for b0, d in sg.hierarchy_is_ordinary().items():
        lines.append(f"    b0={b0:4.1f}:  g_T^2={d['g_T^2']:.4f}   alpha_T={d['alpha_T']:.4f}")
    mat, _, _ = sg.rescaled_matrix_GeV()
    lines += ["", "  RESCALED 3x4 matrix (transmuted M_Pl -> Lambda), GeV:"]
    for t in sp.TOWERS:
        lines.append(f"    {t:16s}" + "".join(f"{v:12.3e}" for v in mat[t]))
    lines += [
        "",
        f"  rescaled top = {s['rescaled_top_GeV']:.0f} GeV (obs 173); the heavy end is now in",
        f"  the right ballpark. The electron at {s['rescaled_electron_GeV']:.1f} GeV (obs 0.5e-3)",
        f"  is the STRUCTURE gap (generation spread), now cleanly separated from the scale.",
        "",
        "READING: the 10^17 scale gap is dimensional transmutation -- exp(-2pi/(b0 g_T^2)) --",
        "of a propagating spin-gravity coupling g_T, the 'other coupling' beyond Newton's G.",
        "An ORDINARY g_T^2~O(0.01-0.1), asymptotically free, condenses at the electroweak scale,",
        "exactly as Lambda_QCD sits 10^19 below M_Pl from an O(1/30) coupling. The torsiton",
        "binds at Lambda, not M_Pl; the Weltformel matrix rescales into the observed ballpark.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "scale_gap.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    # ---- A: the transmutation curve + running ----
    g2 = np.linspace(0.005, 0.5, 400)
    ratio = np.exp(-2 * np.pi / (7.0 * g2))
    axA.semilogy(g2, ratio, "C0", lw=2.2, label=r"$\Lambda/M_{Pl}=e^{-2\pi/(b_0 g_T^2)}$ ($b_0$=7)")
    target = Lam / M_Pl
    g2_hit = sg.g2_for_scale(b0=7.0)
    axA.axhline(target, color="C3", ls="--", lw=1.2,
                label=f"electroweak/Planck ~ {target:.0e}")
    axA.axvline(g2_hit, color="C3", ls=":", lw=1.0)
    axA.scatter([g2_hit], [target], color="C3", s=80, zorder=5, edgecolor="k")
    axA.annotate(f"ordinary coupling\n$g_T^2={g2_hit:.3f}$\n($\\alpha_T={g2_hit/(4*np.pi):.4f}$)",
                 (g2_hit, target), textcoords="offset points", xytext=(35, 20), fontsize=9,
                 color="C3", arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    axA.set_xlabel(r"spin-gravity coupling $g_T^2$")
    axA.set_ylabel(r"$\Lambda / M_{Pl}$ (transmuted scale)")
    axA.set_ylim(1e-20, 2)
    axA.set_title("The scale gap is $\\exp(-40)$, not a tuning\n"
                  "an ordinary $g_T^2$ transmutes Planck $\\to$ electroweak", fontsize=11)
    axA.legend(fontsize=8.5, loc="lower right")
    axA.grid(True, which="both", alpha=0.15)

    # ---- B: matrix before/after transmutation ----
    fm = sp.fermion_matrix_GeV()
    for t in sp.TOWERS:
        before = fm["matrix_GeV"][t]
        after = mat[t]
        obs = sp.OBSERVED_GEV[t]
        axB.scatter([0, 1, 2], before, color=TCOLOR[t], marker="^", s=45, alpha=0.5)
        axB.scatter([0, 1, 2], after, color=TCOLOR[t], marker="s", s=55, edgecolor="k", lw=0.5)
        axB.scatter([0, 1, 2], obs, color=TCOLOR[t], marker="o", s=70, edgecolor="k", lw=0.6,
                    label=t)
    axB.axhspan(1e15, 1e20, color="0.85", alpha=0.5)
    axB.text(0.05, 3e18, "Planck (before): triangles", fontsize=8.5, color="0.4")
    axB.text(0.05, 4e2, "transmuted: squares   observed: circles", fontsize=8.5, color="0.3")
    axB.set_yscale("log")
    axB.set_xticks([0, 1, 2]); axB.set_xticklabels(["gen I", "gen II", "gen III"])
    axB.set_ylabel("mass (GeV)")
    axB.set_ylim(1e-12, 1e20)
    axB.set_title("Weltformel matrix: Planck $\\to$ transmuted\n"
                  "scale gap closes (heavy end into the ballpark)", fontsize=11)
    axB.legend(fontsize=7.5, loc="center right")
    axB.grid(True, which="both", axis="y", alpha=0.12)

    fig.suptitle("Closing the scale gap: dimensional transmutation of the propagating "
                 "spin-gravity coupling ($g_T$, the second coupling beyond $G$)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(PDF_DIR, "scale_gap.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "scale_gap.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'scale_gap.pdf')}")


if __name__ == "__main__":
    main()
