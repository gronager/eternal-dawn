#!/usr/bin/env python3
r"""The inhomogeneous shear channel of the extruder.

Anisotropic collapse carries a shear sigma whose density redshifts as a stiff
fluid (rho_sigma ~ a^-6), growing faster than radiation and dominating the
approach to the bounce (BKL/Mixmaster). In plain GR the shear-viscous entropy
production diverges at the singularity. The Einstein--Cartan torsion cap bounds
the shear at the bounce scale, so the entropy integral J_shear stays O(1) and

    ln D_shear = 2 eta~ (Omega_bounce/T_bounce) J_shear

is suppressed by the same ~1e-11 prefactor as the bulk channel: the bounce is
adiabatic even when shear-dominated. Renders figures/pdf/shear.pdf and writes
sims/output/shear.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import cve_filter as cve
from cartasis_sims import extruder as ex

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    rho_C = 1.0e50
    zt_KSS = 1.0 / (4.0 * np.pi)
    Om = cve.bounce_rate_GeV(rho_C)
    T_b = cve.bounce_temperature_GeV(rho_C)

    fracs = [0.0, 0.1, 0.5, 0.9, 0.99]
    rows = []
    for fs in fracs:
        J = ex.shear_entropy_integral(fs)
        Dk = ex.shear_dilution_factor(zt_KSS, fs) - 1.0
        rows.append((fs, J, Dk))

    lines = [
        "The inhomogeneous shear channel of the extruder",
        "=" * 58,
        f"  bounce rate Omega = {Om:.2e} GeV,  T_bounce = {T_b:.2e} GeV",
        f"  prefactor Omega/T = {Om/T_b:.2e}   (shared with the bulk channel)",
        "",
        "  shear-dominated bounces stay adiabatic (eta~ = 1/4pi, KSS):",
        "    f_sigma   J_shear     D - 1",
    ]
    for fs, J, Dk in rows:
        lines.append(f"    {fs:5.2f}    {J:7.4f}    {Dk:.3e}")
    lines += [
        "",
        "READING: the shear density rho_sigma ~ a^-6 outgrows radiation and would,",
        "in GR, drive sigma^2/T -> infinity at the singularity -- a divergent",
        "(BKL) entropy production. The torsion bounce caps the total density at",
        "rho_C, so the shear is bounced too: sigma^2 <= rho_C and J_shear stays",
        "O(1) even for a maximally shear-dominated bounce (f_sigma -> 1). With the",
        "shared 1e-11 prefactor, D - 1 < 1e-10 across the whole range. The same",
        "torsion that removes the singularity tames the BKL entropy divergence:",
        "inheritance survives the inhomogeneous channel too.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "shear.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # Panel L: a shear-dominated bounce -- shear is capped, not divergent.
    tau, a, H, rho_sig = ex.simulate_anisotropic_bounce(f_sigma=0.9)
    # mirror to the collapse branch for a symmetric picture
    tau_full = np.concatenate([-tau[::-1], tau])
    a_full = np.concatenate([a[::-1], a])
    sig2_full = np.concatenate([rho_sig[::-1], rho_sig])
    kernel = sig2_full * a_full
    m = np.abs(tau_full) < 5.0
    axL.plot(tau_full[m], a_full[m], "C0", lw=1.8, label=r"scale factor $a$")
    axL.plot(tau_full[m], sig2_full[m], "C3", lw=1.8,
             label=r"shear $\sigma^2\propto a^{-6}$ (capped at $\rho_C$)")
    axL.fill_between(tau_full[m], 0, kernel[m], color="C2", alpha=0.35,
                     label=r"entropy kernel $\sigma^2 a$")
    axL.axhline(1.0, color="0.6", ls=":", lw=1.0)
    axL.text(-4.8, 1.03, r"$\rho_C$ cap", fontsize=8, color="0.4")
    axL.annotate("GR: $\\sigma^2\\to\\infty$\n(BKL, no bounce)", xy=(0, 0.9),
                 xytext=(1.4, 2.4), fontsize=8, color="0.35",
                 arrowprops=dict(arrowstyle="->", color="0.5", lw=1.0))
    axL.set_xlabel(r"dimensionless proper time  $\tau\,\sqrt{8\pi G\rho_C/3}$")
    axL.set_ylabel("dimensionless")
    axL.set_title("A shear-dominated bounce ($f_\\sigma=0.9$):\n"
                  "torsion caps the shear instead of letting it diverge",
                  fontsize=11)
    axL.set_ylim(0, 3.2)
    axL.legend(fontsize=8, loc="upper right")
    axL.grid(True, alpha=0.2)

    # Panel R: dilution vs anisotropy fraction -- adiabatic across the board.
    fs_grid = np.linspace(0.0, 0.99, 60)
    for zt, lab, col in [(zt_KSS, r"$\tilde\eta=1/4\pi$ (KSS)", "C0"),
                         (1.0, r"$\tilde\eta=1$", "C1")]:
        Dm1 = [ex.shear_dilution_factor(zt, fs) - 1.0 for fs in fs_grid]
        axR.plot(fs_grid, Dm1, col, lw=2.0, label=lab)
    D_hor = cve.dilution_horizon(1.0e53, rho_C) - 1.0
    axR.axhline(D_hor, color="0.4", ls="--", lw=1.2)
    axR.text(0.02, D_hor * 0.03, r"horizon-saturating $D_{\max}\sim10^{49}$",
             fontsize=8, color="0.3")
    axR.set_yscale("log")
    axR.set_xlabel(r"shear fraction at the bounce  $f_\sigma$")
    axR.set_ylabel(r"entropy dilution  $D-1$")
    axR.set_title("Even shear-dominated, the bounce stays adiabatic\n"
                  r"$D-1\lesssim10^{-10}\Rightarrow\eta_{\rm baby}=\eta_{\rm parent}$",
                  fontsize=11)
    axR.set_ylim(1e-13, 1e52)
    axR.legend(fontsize=9, loc="center right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "shear.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
