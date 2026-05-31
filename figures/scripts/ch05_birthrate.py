#!/usr/bin/env python3
r"""Pinning the OG birth rate beta as a de Sitter nucleation rate.

beta = lambda * H_Lambda^4/c^3: the PREFACTOR is fixed by the observed dark energy
(~4e-97 m^-3 s^-1), the dimensionless lambda = exp(-I) with I the bounce instanton
action. The percolation threshold lambda_crit ~ 0.24 (I_crit ~ 1.4) splits the
supraverse: packed polygonal foam (I < 1.4) vs dilute island universes (I > 1.4).
Renders figures/pdf/birth_rate.pdf, writes sims/output/birth_rate.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import birth_rate as br

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    pref = br.nucleation_prefactor()
    H_L = br.hubble_lambda()
    I_crit = br.instanton_action(br.LAMBDA_CRIT)

    lines = [
        "Pinning the OG birth rate beta (as a de Sitter nucleation rate)",
        "=" * 60,
        f"  H_Lambda = H_0 sqrt(Omega_L) = {H_L:.3e} s^-1   (the void's dark energy)",
        f"  prefactor H_Lambda^4/c^3     = {pref:.3e} m^-3 s^-1   <- PINNED by Lambda",
        f"  beta = lambda * prefactor,   lambda = exp(-I),  I = bounce instanton action",
        "",
        f"  percolation threshold: lambda_crit = {br.LAMBDA_CRIT}, "
        f"I_crit = {I_crit:.2f}",
        "    lambda   I=-ln(lam)   beta[/m^3/s]    sep[R_H]    regime",
    ]
    for lam in (1.0, br.LAMBDA_CRIT, 1e-3, 1e-20):
        lines.append(f"    {lam:.0e}    {br.instanton_action(lam):5.1f}     "
                     f"{br.beta(lam):.2e}   {br.ogu_separation_hubble(lam):.1e}    "
                     f"{'packed' if br.percolates(lam) else 'dilute'}")
    lines += [
        "",
        "READING: beta is no longer a free parameter in SCALE -- its prefactor is the",
        "de Sitter attempt rate set by the observed dark energy (~4e-97 m^-3 s^-1, one",
        "per Hubble 4-volume). The one remaining unknown is the dimensionless lambda =",
        "exp(-I), I the rho_C-crossing gravitational instanton action (Q4). But lambda",
        "has a SHARP critical value: de Sitter percolation at lambda_crit ~ 0.24",
        "(I_crit ~ 1.4) decides the supraverse's character --",
        "  I < 1.4 (near-Planckian seed): OGUs percolate -> PACKED polygonal foam that",
        "           coarsens (the void_foam / coarsening picture);",
        "  I > 1.4 (macroscopic seed): the void inflates between births -> DILUTE,",
        "           isolated, round island universes (eternal-inflation-like, no",
        "           coarsening).",
        "We are a BHU1 inside an OGU either way (the internal black-hole channel is",
        "independent of OGU spacing). So: prefactor pinned by dark energy, decisive",
        "structure pinned at I_crit~1.4; the regime awaits the bounce instanton.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "birth_rate.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: beta vs lambda, prefactor pinned by dark energy.
    lam = np.logspace(-50, 0, 400)
    axL.plot(lam, br.beta(lam), "C0", lw=2.2)
    axL.axvline(br.LAMBDA_CRIT, color="C3", ls="--", lw=1.3)
    axL.axvspan(br.LAMBDA_CRIT, 1.0, color="0.5", alpha=0.12)
    axL.plot([1.0], [pref], "ko", ms=6)
    axL.annotate(f"prefactor $H_\\Lambda^4/c^3$\n$={pref:.1e}$ m$^{{-3}}$s$^{{-1}}$\n"
                 "(pinned by dark energy)", xy=(1.0, pref), xytext=(1e-30, 1e-100),
                 fontsize=8, arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.text(0.32, 1e-140, "packed\n($\\lambda>\\lambda_{\\rm crit}$)", fontsize=8,
             color="C3", ha="center")
    axL.text(1e-25, 1e-140, "dilute (island universes)", fontsize=8, color="0.4")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"dimensionless nucleation rate  $\lambda=e^{-I}$")
    axL.set_ylabel(r"OG birth rate  $\beta=\lambda\,H_\Lambda^4/c^3$  [m$^{-3}$s$^{-1}$]")
    axL.set_title(r"$\beta$'s scale is fixed by dark energy"
                  "\n"
                  r"only the instanton factor $\lambda$ is free", fontsize=11)
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: regime vs instanton action I = -ln lambda.
    I = np.linspace(0.0, 100.0, 400)
    lam_I = np.exp(-I)
    axR.plot(I, br.ogu_separation_hubble(lam_I), "C0", lw=2.2)
    axR.axvline(I_crit, color="C3", ls="--", lw=1.3)
    axR.axvspan(0.0, I_crit, color="0.5", alpha=0.12)
    axR.text(I_crit + 1.5, 3.0, "$I_{\\rm crit}\\approx1.4$", fontsize=9, color="C3")
    axR.text(0.3, 1e8, "packed\nfoam", fontsize=9, color="C3")
    axR.text(40, 1e3, "dilute island universes\n(eternal-inflation-like)", fontsize=9,
             color="0.4")
    axR.set_yscale("log")
    axR.set_xlabel(r"bounce instanton action  $I=-\ln\lambda$  [$\hbar$]")
    axR.set_ylabel(r"mean OGU separation  [de Sitter horizon radii]")
    axR.set_title("One number sets the supraverse's character\n"
                  "near-Planckian seed packs it; macroscopic seed dilutes it",
                  fontsize=11)
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "birth_rate.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
