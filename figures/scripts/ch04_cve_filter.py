#!/usr/bin/env python3
r"""Does chiral-vortical transport at the rotating bounce source eta ~ 6e-10?

Order-of-magnitude estimate: eta ~ C * (omega/T), with omega/T the
vorticity-to-temperature ratio at the bounce. Renders figures/pdf/cve_filter.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import cve_filter as cve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    rho_C = 1.0e50
    T = cve.bounce_temperature_GeV(rho_C)
    woT = cve.vorticity_over_T(rho_C, spin_fraction=1.0)
    need = cve.spinC_needed(rho_C)

    lines = [
        "Chiral-vortical asymmetry estimate at the rotating bounce",
        "=" * 58,
        f"  Cartasis density rho_C      = {rho_C:.1e} kg/m^3",
        f"  bounce plasma temperature T = {T:.2e} GeV  (EW restored; sphalerons active)",
        f"  bounce rate sqrt(8piG rho/3)= {cve.bounce_rate_s(rho_C):.2e} s^-1",
        f"  omega/T (maximal spin)      = {woT:.2e}",
        f"  eta ~ C * (omega/T), observed eta_obs = {cve.ETA_OBS:.1e}",
        f"  => need  C * spin_fraction  = {need:.1f} to match eta_obs",
        "",
        "READING: omega/T ~ 1e-11 at the bounce sits within ~1-2 orders of the",
        "observed eta ~ 6e-10. The mechanism is the right ballpark (it could have",
        "been 1e-30 or 1e+5); it generically UNDERSHOOTS by ~1-2 orders, so it",
        "needs near-maximal bounce vorticity and an O(10) anomaly/sphaleron factor.",
        "Robust qualitative result: the asymmetry is sourced ALONG the spin axis",
        "and scales with bounce vorticity -- tying baryogenesis to the Kerr/axis story.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "cve_filter.txt"), "w") as f:
        f.write(text + "\n")

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    fs = np.linspace(0.0, 1.0, 200)
    for C, col, ls in [(0.1, "C0", "--"), (1.0, "C0", "-"), (10.0, "C0", "-.")]:
        eta = [cve.asymmetry_estimate(rho_C, spin_fraction=f, C=C) for f in fs]
        ax.plot(fs, eta, col, ls=ls, lw=1.8, label=fr"$C={C:g}$")
    ax.axhspan(3e-10, 9e-10, color="C3", alpha=0.18)
    ax.axhline(cve.ETA_OBS, color="C3", lw=1.4)
    ax.text(0.02, cve.ETA_OBS * 1.25, r"observed $\eta\simeq6\times10^{-10}$",
            color="C3", fontsize=9)
    ax.set_yscale("log")
    ax.set_xlabel("bounce spin fraction  $f_\\omega=\\omega/\\Omega_{\\rm bounce}$")
    ax.set_ylabel(r"chiral asymmetry estimate  $\eta\sim C\,\omega/T$")
    ax.set_title("Chiral-vortical asymmetry at the rotating bounce\n"
                 fr"$T\simeq{T:.0e}$ GeV, $\omega/T|_{{\max}}\simeq{woT:.0e}$",
                 fontsize=11)
    ax.set_ylim(1e-14, 1e-8)
    ax.legend(title="anomaly/sphaleron factor", fontsize=9)
    ax.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "cve_filter.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
