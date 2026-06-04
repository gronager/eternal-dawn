#!/usr/bin/env python3
r"""The bounce coupling is gravitational, not the mass coupling -- a neutron star decides.

The Einstein--Cartan bounce density scales inversely with the four-fermion coupling. Plotted on
a density axis: the gravitational-torsion bounce (Part I, ~1e50 kg/m^3); where the strong,
mass-giving coupling (~1/Lambda^2, ~10^33 stronger) would instead put it (~0.2x nuclear); and the
neutron-star core band that already sits ABOVE that without bouncing -- excluding the strong
coupling as the bounce driver. Hence the bounce and the mass scale are two separate couplings, and
Part I's bounce is undisturbed by the mass sector (melted to radiation at the bounce temperature).

Renders figures/pdf/bounce_consistency.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import bounce_consistency as bc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
for d in (PDF_DIR, OUT_DIR):
    os.makedirs(d, exist_ok=True)


def main():
    s = bc.summary()
    nst = s["neutron_star_test"]
    lines = [
        "The bounce coupling is gravitational, not the mass coupling",
        "=" * 58,
        f"  mass coupling / gravitational coupling = (M_Pl/Lambda)^2 = "
        f"{s['coupling_ratio_strong_over_grav']:.2e}",
        f"  bounce density (gravitational, Part I): {bc.RHO_C_GRAV:.0e} kg/m^3",
        f"  bounce density IF driven by the strong/mass coupling: "
        f"{nst['rho_C_strong_kg_m3']:.2e} kg/m^3 = {nst['in_nuclear_units']:.2f}x nuclear",
        f"  neutron-star core ~ {bc.RHO_NS_CORE:.0e} kg/m^3 exceeds that by "
        f"{nst['ns_core_exceeds_by']:.0f}x WITHOUT bouncing -> strong coupling EXCLUDED",
        f"  bounce temperature ~ {s['bounce_temperature_GeV']:.1e} GeV; "
        f"T_bounce/T_c ~ {s['T_bounce_over_T_c']:.0e} (condensate melted at the bounce)",
        "",
        "READING: " + s["verdict"],
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "bounce_consistency.txt"), "w") as f:
        f.write(text + "\n")

    fig, ax = plt.subplots(figsize=(11.0, 4.4))
    # density axis (log), kg/m^3
    ax.set_xscale("log")
    ax.set_xlim(1e14, 1e52)
    ax.set_ylim(0, 1)
    ax.get_yaxis().set_visible(False)

    # reference bands
    ax.axvspan(bc.RHO_NUCLEAR, bc.RHO_NS_CORE, color="C0", alpha=0.18)
    ax.axvline(bc.RHO_NUCLEAR, color="C0", ls=":", lw=1.0)
    ax.text(bc.RHO_NUCLEAR, 0.86, "  nuclear density", fontsize=8.5, color="C0")
    ax.text(np.sqrt(bc.RHO_NUCLEAR * bc.RHO_NS_CORE), 0.07, "neutron-star cores\n(do NOT bounce)",
            fontsize=8.5, color="C0", ha="center")

    # the strong-coupling bounce (excluded)
    rho_s = nst["rho_C_strong_kg_m3"]
    ax.scatter([rho_s], [0.5], s=160, marker="X", color="C3", zorder=5, edgecolor="k")
    ax.annotate("if the MASS coupling drove the bounce:\n"
                f"$\\rho_C\\sim{rho_s:.0e}$ kg/m$^3$ (0.2$\\times$ nuclear)\n"
                "-- excluded: neutron stars sit above it",
                (rho_s, 0.5), textcoords="offset points", xytext=(10, 28), fontsize=8.5,
                color="C3", arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))

    # the gravitational bounce (Part I)
    ax.scatter([bc.RHO_C_GRAV], [0.5], s=170, marker="o", color="C2", zorder=5, edgecolor="k")
    ax.annotate("the actual bounce (gravitational torsion,\n"
                f"$G/c^4$): $\\rho_C\\sim10^{{50}}$ kg/m$^3$",
                (bc.RHO_C_GRAV, 0.5), textcoords="offset points", xytext=(-30, -42),
                fontsize=8.5, color="C2", ha="center",
                arrowprops=dict(arrowstyle="->", color="C2", lw=0.8))

    # arrow showing the factor between them
    ax.annotate("", xy=(bc.RHO_C_GRAV, 0.5), xytext=(rho_s, 0.5),
                arrowprops=dict(arrowstyle="<->", color="0.5", lw=1.2, ls="--"))
    ax.text(np.sqrt(rho_s * bc.RHO_C_GRAV), 0.56,
            f"$\\times{s['coupling_ratio_strong_over_grav']:.0e}$\n(coupling ratio)",
            fontsize=8, color="0.4", ha="center")

    ax.set_xlabel(r"density at which the collapse bounces (kg m$^{-3}$)")
    ax.set_title("Two couplings, and a neutron star proves it: the bounce is gravitational, "
                 "the mass scale a separate $\\sim10^{33}$ stronger sector", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "bounce_consistency.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "bounce_consistency.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'bounce_consistency.pdf')}")


if __name__ == "__main__":
    main()
