#!/usr/bin/env python3
r"""Computing the OG-nucleation instanton action I -- and the verdict it gives.

I = 2 pi M_seed c^2/(hbar H_Lambda): the de Sitter Boltzmann / horizon-deficit action
to spike the cold void to rho_C. For the minimal bounce seed (~9e14 kg) and our dark
energy, I ~ 2.5e84 -- 84 orders above the percolation threshold I_crit ~ 1.4. So the
supraverse is DILUTE: isolated, round, eternal-inflation-like island universes, not a
packed foam. Renders figures/pdf/instanton.pdf, writes sims/output/instanton.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import birth_rate as br
from cartasis_sims import cycle as cy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

M_PLANCK = 2.176e-8       # kg


def main() -> None:
    M_seed = br.minimal_seed_mass()
    I = br.seed_instanton_action()
    Icrit = br.instanton_action_threshold()
    T_dS = cy.de_sitter_temperature(br.hubble_lambda())
    M_thresh = br.seed_mass_for_threshold()
    T_for_packed = M_seed * 3e8**2 / (1.380649e-23 * Icrit)

    lines = [
        "Computing the OG-nucleation instanton action I",
        "=" * 58,
        "  I = 2 pi M_seed c^2/(hbar H_Lambda) = M_seed c^2/(k_B T_dS)",
        "    (de Sitter Boltzmann factor = cosmological-horizon entropy deficit)",
        "",
        f"  minimal bounce seed M_seed   = {M_seed:.2e} kg  (R_s = R at rho_C, ~1e-12 m)",
        f"  de Sitter void temperature   = {T_dS:.2e} K",
        f"  => instanton action I        = {I:.2e}",
        f"  percolation threshold I_crit = {Icrit:.2f}",
        f"  I / I_crit                   = {I/Icrit:.1e}",
        "",
        f"  to get I = I_crit you would need:",
        f"    a seed of {M_thresh:.1e} kg  (sub-Planckian by ~{np.log10(M_PLANCK/M_thresh):.0f}"
        f" orders -- cannot self-gravitate or bounce), OR",
        f"    a void at {T_for_packed:.1e} K  (vs the actual {T_dS:.0e} K -- 84 orders too cold).",
        "",
        "VERDICT: I ~ 1e84 >> I_crit ~ 1.4, robust by 84 orders of magnitude (no",
        "semiclassical correction flips it). So lambda = exp(-I) ~ 0 and the supraverse",
        "is DILUTE: OGUs are isolated, round, eternal-inflation-like island universes,",
        "astronomically far apart -- NOT a packed polygonal foam.",
        "",
        "Reconciliation: the packed-foam constructions (average size by impingement,",
        "radiative coarsening, the polygonal render) were the PERCOLATING limit; the",
        "physical regime is dilute. What survives untouched: we are a BHU1 inside an",
        "OGU (the internal black-hole channel is independent of OGU spacing); the",
        "recursive cycle (each isolated island heat-dies and seeds its own dilute",
        "scatter); M_OGU ~ c^3 t/2G; and the dark-energy-pinned prefactor. The dilute",
        "answer also vindicates the earlier reading that well-separated universes are",
        "round, not polygonal.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "instanton.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: I vs seed mass.
    M = np.logspace(-72, 20, 400)
    axL.plot(M, br.de_sitter_boltzmann_action(M), "C0", lw=2.2)
    axL.axhline(Icrit, color="C3", ls="--", lw=1.3)
    axL.text(1e-70, 4.0, r"$I_{\rm crit}\approx1.4$ (percolation)", fontsize=8,
             color="C3")
    axL.axvline(M_PLANCK, color="0.5", ls=":", lw=1.0)
    axL.text(M_PLANCK * 2, 1e60, "Planck mass", fontsize=8, color="0.4", rotation=90)
    axL.plot([M_seed], [I], "ko", ms=6)
    axL.annotate(f"minimal bounce seed\n$M_{{\\rm seed}}\\sim10^{{15}}$ kg\n"
                 f"$I\\sim10^{{84}}$", xy=(M_seed, I), xytext=(1e-30, 1e70),
                 fontsize=8, arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.axhspan(Icrit, 1e95, color="0.5", alpha=0.10)
    axL.text(1e-60, 1e30, "DILUTE\n(island universes)", fontsize=9, color="0.35")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"seed mass  $M_{\rm seed}$  [kg]")
    axL.set_ylabel(r"instanton action  $I=2\pi M_{\rm seed}c^2/\hbar H_\Lambda$")
    axL.set_title("Any seed that can bounce gives $I\\gg I_{\\rm crit}$\n"
                  "(the $I_{\\rm crit}$ seed is sub-Planckian and cannot bounce)",
                  fontsize=11)
    axL.set_ylim(1e-2, 1e95)
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: I vs void temperature.
    T = np.logspace(-32, 58, 400)
    I_of_T = M_seed * 3e8**2 / (1.380649e-23 * T)
    axR.plot(T, I_of_T, "C0", lw=2.2)
    axR.axhline(Icrit, color="C3", ls="--", lw=1.3)
    axR.axvline(T_dS, color="C0", ls=":", lw=1.2)
    axR.plot([T_dS], [I], "ko", ms=6)
    axR.annotate(f"our cold void\n$T_{{\\rm dS}}\\sim10^{{-30}}$ K\n$I\\sim10^{{84}}$",
                 xy=(T_dS, I), xytext=(1e-22, 1e40), fontsize=8,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axR.axvline(T_for_packed, color="C3", ls=":", lw=1.0)
    axR.text(T_for_packed * 0.5, 1e6, "would need\n$\\sim10^{54}$ K\nfor packed",
             fontsize=8, color="C3", ha="right")
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"void temperature  $T_{\rm dS}$  [K]")
    axR.set_ylabel(r"instanton action  $I=M_{\rm seed}c^2/k_B T_{\rm dS}$")
    axR.set_title("The void is 84 orders too cold to percolate\n"
                  "$\\Rightarrow$ dilute supraverse", fontsize=11)
    axR.set_ylim(1e-2, 1e95)
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "instanton.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
