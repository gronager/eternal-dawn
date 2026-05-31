#!/usr/bin/env python3
r"""The extruder: the entropy budget of the Einstein--Cartan bounce.

    "Universes are not born, they are extruded."

How much entropy does the bounce generate? The dilution factor D = S_after/S_before
sets eta_baby = eta_parent / D and decides whether a baby inherits its parent's
asymmetry (D~1, pure lineages) or has it crushed (D large). We compute D from the
bounce dynamics with a bulk viscosity zeta = zeta_tilde * s:

    ln D = 9 zeta_tilde (Omega_bounce/T_bounce) J,   J = int H_dim^2 a dtau_dim.

Renders figures/pdf/extruder.pdf and writes sims/output/extruder.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import bounce as bnc
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
    res = ex.analyze(w=1.0 / 3.0, rho_C=rho_C)

    # viscosity needed to reach the horizon-saturating dilution (absurd)
    D_hor = cve.dilution_horizon(1.0e53, rho_C)
    zt_hor = np.log(D_hor) / (9.0 * res.adiabaticity * res.J)
    zt_KSS = 1.0 / (4.0 * np.pi)        # a shear-viscosity reference scale

    lines = [
        "The extruder: entropy budget of the Einstein-Cartan bounce",
        "=" * 58,
        '  "Universes are not born, they are extruded."',
        "",
        f"  Cartasis density rho_C        = {rho_C:.1e} kg/m^3   (radiation, w=1/3)",
        f"  bounce temperature T_bounce   = {res.T_bounce_GeV:.2e} GeV",
        f"  bounce rate Omega_bounce      = {res.Omega_bounce_GeV:.2e} GeV",
        f"  adiabaticity Omega/T          = {res.adiabaticity:.2e}   (= omega/T of the CVE)",
        f"  entropy-production integral J = {res.J:.3f}   (dimensionless, O(1))",
        "",
        "  ln D = 9 zeta~ (Omega/T) J,  for bulk viscosity zeta = zeta~ * s:",
        f"    zeta~ = 0      (conformal radiation) -> D = 1 exactly  (zeta=0)",
        f"    zeta~ = 1/4pi  (KSS shear scale)     -> D - 1 = {ex.dilution_factor(zt_KSS, res)-1:.2e}",
        f"    zeta~ = 1      (O(1) dissipation)    -> D - 1 = {ex.dilution_factor(1.0, res)-1:.2e}",
        f"    zeta~ = 1e3    (strong dissipation)  -> D - 1 = {ex.dilution_factor(1e3, res)-1:.2e}",
        "",
        f"  To reach the horizon-saturating D = {D_hor:.1e} would require",
        f"    zeta~ = {zt_hor:.1e}  -- unphysical by ~12 orders of magnitude.",
        "",
        "READING: the prefactor Omega/T ~ 1e-11 -- the SAME small parameter that",
        "limits the chiral-vortical estimate -- means the bounce is too FAST to",
        "thermalize: it is adiabatic to ~1 part in 1e11 even for O(1) viscosity,",
        "and a conformal radiation fluid generates no bulk-viscous entropy at all.",
        "So D = 1 to extraordinary precision: eta_baby = eta_parent. Inheritance",
        "is robust and the horizon-saturating limit is dynamically out of reach.",
        "The factor-of-41 is the OGU's bill; every descendant simply inherits the",
        "~1 ppb its parent already made, and the founder's chiral sign is protected",
        "-> PURE lineages. (Homogeneous channel; inhomogeneous shear / particle",
        "production at the bounce remain to be checked -- Tier 1.)",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "extruder.txt"), "w") as f:
        f.write(text + "\n")

    sol = bnc.simulate_bounce(w=1.0 / 3.0, a_init=6.0, t_max=30.0, n_points=8000)
    kernel = sol.H**2 * sol.a       # entropy-production integrand (area = J)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # Panel L: where (little) entropy is generated through the bounce.
    m = np.abs(sol.t) < 6.0
    axL.plot(sol.t[m], sol.a[m], "C0", lw=1.8, label=r"scale factor $a$")
    axL.plot(sol.t[m], sol.H[m], "C1", lw=1.4, label=r"Hubble rate $H$")
    axL.fill_between(sol.t[m], 0, kernel[m], color="C2", alpha=0.35,
                     label=r"entropy kernel $H^2 a$")
    axL.plot(sol.t[m], kernel[m], "C2", lw=1.6)
    axL.axvline(0.0, color="0.6", ls=":", lw=1.0)
    axL.text(0.15, 1.05, "bounce", fontsize=8, color="0.4")
    axL.set_xlabel(r"dimensionless proper time  $\tau\,\sqrt{8\pi G\rho_C/3}$")
    axL.set_ylabel("dimensionless")
    axL.set_title(r"Entropy is generated on the shoulders, not at $H{=}0$"
                  "\n"
                  fr"$\int H^2 a\,d\tau = J = {res.J:.2f}$  (area shaded)",
                  fontsize=11)
    axL.legend(fontsize=9, loc="upper right")
    axL.grid(True, alpha=0.2)

    # Panel R: dilution vs viscosity -- adiabatic plateau, horizon target far off.
    zts = np.logspace(-3.0, 12.0, 400)
    Dm1 = np.array([ex.dilution_factor(z, res) - 1.0 for z in zts])
    axR.plot(zts, np.clip(Dm1, 1e-20, None), "C3", lw=2.0,
             label=r"$D-1$")
    axR.axvspan(1e-3, 1.0, color="C0", alpha=0.15)
    axR.text(3e-3, 1e-2, "physical\nviscosity\n$\\tilde\\zeta\\lesssim1$",
             fontsize=8, color="C0")
    axR.axhline(D_hor - 1.0, color="0.4", ls="--", lw=1.2)
    axR.text(1e-3, (D_hor - 1.0) * 2.0,
             r"horizon-saturating $D_{\max}\sim10^{49}$", fontsize=8, color="0.3")
    axR.axvline(zt_hor, color="0.4", ls=":", lw=1.0)
    axR.text(zt_hor * 0.18, 1e-15,
             fr"need $\tilde\zeta\sim10^{{{int(round(np.log10(zt_hor)))}}}$",
             fontsize=8, color="0.3", ha="right", rotation=90)
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"bulk viscosity  $\tilde\zeta=\zeta/s$")
    axR.set_ylabel(r"entropy dilution  $D-1$")
    axR.set_title("The bounce is adiabatic for any physical viscosity\n"
                  r"$\eta_{\rm baby}=\eta_{\rm parent}/D\Rightarrow$ inheritance",
                  fontsize=11)
    axR.set_ylim(1e-18, 1e52)
    axR.legend(fontsize=9, loc="lower right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "extruder.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
