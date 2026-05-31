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

    # Inheritance / dilution track: our parent (~10^53 kg) already sits at
    # eta_parent ~ eta_obs; the baby inherits eta_parent / D.
    M_parent = 1.0e53
    D_hor = cve.dilution_horizon(M_parent, rho_C)
    # The dilution D that reproduces the observed eta from an eta_obs parent:
    D_obs = cve.ETA_OBS / cve.ETA_OBS  # = 1 by construction (parent==obs anchor)

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
        "TRACK A -- OGU must MANUFACTURE eta (no parent to inherit from):",
        "  omega/T ~ 1e-11 at the bounce sits within ~1-2 orders of the observed",
        "  eta ~ 6e-10. Right ballpark (could have been 1e-30 or 1e+5); generically",
        "  UNDERSHOOTS by ~41x, so the OGU needs near-maximal bounce vorticity and",
        "  an O(10) anomaly/sphaleron factor C. Sourced ALONG the spin axis, tying",
        "  baryogenesis to the Kerr/axis story.",
        "",
        "TRACK B -- BHU_n INHERITS eta (our case, by the DM+DE anchor):",
        "  sphalerons conserve B-L, so the parent's net B-L passes through as a",
        "  protected remnant. eta_baby = eta_parent / D, with D>=1 the bounce",
        "  entropy-dilution factor. Parent already sits at eta_parent ~ eta_obs, so:",
        f"    D = 1     (adiabatic extrusion)      -> eta_baby = {cve.inherited_eta(dilution=1.0):.1e}  (= eta_obs)",
        f"    D = 1e2   (mild irreversible heating) -> eta_baby = {cve.inherited_eta(dilution=1e2):.1e}",
        f"    D = {D_hor:.1e} (horizon-saturating)  -> eta_baby = {cve.inherited_eta(dilution=D_hor):.1e}",
        f"  [horizon limit D_max = S_BH/S_rad for M_parent={M_parent:.0e} kg, T at bounce]",
        "",
        "READING: the factor-of-41 is the OGU's bill, not ours -- a BHU inherits",
        "the ~1 ppb its parent already made. The crux is the dilution D. The two",
        "anchors are D~1 (adiabatic) and D~1e49 (horizon-saturating); the latter",
        "would crush eta to ~1e-59, which we DO NOT observe. So the measured",
        "eta~6e-10 itself selects the near-adiabatic end (D within a few orders of",
        "1): the extrusion does NOT thermalize to the horizon entropy. That same",
        "low-D regime means inheritance dominates the fresh skew -> the founder's",
        "sign is protected -> lineages stay chirally PURE. High-D would force fresh",
        "skew, re-roll the sign, and MIX lineages. Pinning D needs the real bounce",
        "entropy budget (open question Q16a / Q8).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "cve_filter.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # Panel A: Track A -- fresh manufacture vs spin (the OGU's problem).
    fs = np.linspace(0.0, 1.0, 200)
    for C, ls in [(0.1, "--"), (1.0, "-"), (10.0, "-.")]:
        eta = [cve.asymmetry_estimate(rho_C, spin_fraction=f, C=C) for f in fs]
        axA.plot(fs, eta, "C0", ls=ls, lw=1.8, label=fr"$C={C:g}$")
    axA.axhspan(3e-10, 9e-10, color="C3", alpha=0.18)
    axA.axhline(cve.ETA_OBS, color="C3", lw=1.4)
    axA.text(0.02, cve.ETA_OBS * 1.25, r"observed $\eta\simeq6\times10^{-10}$",
             color="C3", fontsize=9)
    axA.set_yscale("log")
    axA.set_xlabel(r"bounce spin fraction  $f_\omega=\omega/\Omega_{\rm bounce}$")
    axA.set_ylabel(r"asymmetry  $\eta\sim C\,\omega/T$")
    axA.set_title(f"Track A: OGU manufactures $\\eta$\n"
                  fr"$T\simeq{T:.0e}$ GeV, $\omega/T|_{{\max}}\simeq{woT:.0e}$",
                  fontsize=11)
    axA.set_ylim(1e-14, 1e-8)
    axA.legend(title="anomaly/sphaleron $C$", fontsize=9)
    axA.grid(True, which="both", alpha=0.2)

    # Panel B: Track B -- inherited eta vs entropy-dilution D (our case).
    Ds = np.logspace(0.0, 50.0, 300)
    eta_inh = [cve.inherited_eta(cve.ETA_OBS, D) for D in Ds]
    axB.plot(Ds, eta_inh, "C2", lw=2.0,
             label=r"$\eta_{\rm baby}=\eta_{\rm parent}/D$")
    axB.axhline(cve.asymmetry_estimate(spin_fraction=1.0, C=1.0), color="C0",
                ls=":", lw=1.4, label=r"fresh skew floor ($C{=}1$, max spin)")
    axB.axhspan(3e-10, 9e-10, color="C3", alpha=0.18)
    axB.axhline(cve.ETA_OBS, color="C3", lw=1.4)
    axB.axvline(1.0, color="0.4", ls="--", lw=1.0)
    axB.text(2.5, 1e-30, "adiabatic\n$D\\simeq1$", fontsize=8, color="0.3")
    axB.axvline(D_hor, color="0.4", ls="--", lw=1.0)
    axB.text(D_hor * 0.4, 1e-30, "horizon-\nsaturating\n$D_{\\max}$",
             fontsize=8, color="0.3", ha="right")
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlabel(r"bounce entropy-dilution factor  $D=S_{\rm after}/S_{\rm before}$")
    axB.set_ylabel(r"inherited asymmetry  $\eta_{\rm baby}$")
    axB.set_title("Track B: BHU inherits $\\eta$\n"
                  "observation selects the near-adiabatic end", fontsize=11)
    axB.set_ylim(1e-62, 1e-7)
    axB.legend(fontsize=8, loc="lower left")
    axB.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "cve_filter.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
