#!/usr/bin/env python3
r"""Gravity-scaling magnification: why the wallpaper is two-tier.

Gravity-scale every object to a common reference density; its magnification relative to
true size is (rho/rho_ref)^(1/3). The result is two-tier:
 * an OGU is a horizon-mass hole with ~ the cosmic/void density (R_s = R_H), so its
   magnification vs the void is ~1 -- gravity-scaled coords ~ TRUE coords, the dilute
   scatter survives, the OGU-level frame is honestly near-empty;
 * astrophysical black holes are vastly denser than their host universe, so their
   magnification is ~1e9-1e14 -- the conformal map blows them up so nested structure
   shows. Renders figures/pdf/gravity_scale.pdf, writes sims/output/gravity_scale.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import gravity_scale as gs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    objs = [
        ("OGU (horizon mass)\nvs void", 9.2e52, gs.RHO_LAMBDA, "C0"),
        ("SMBH $10^9 M_\\odot$\nvs universe", 2e39, gs.RHO_CRIT, "0.5"),
        ("SMBH $10^6 M_\\odot$\nvs universe", 2e36, gs.RHO_CRIT, "0.5"),
        ("stellar $10 M_\\odot$\nvs universe", 2e31, gs.RHO_CRIT, "0.5"),
    ]
    lines = [
        "Gravity-scaling magnification: a two-tier wallpaper",
        "=" * 58,
        "  mu = (rho_object / rho_reference)^(1/3); rho_ref cancels in any ratio.",
        "",
        "  object                         rho_mean [kg/m^3]   magnification mu",
    ]
    mags = []
    for label, M, rho_ref, _ in objs:
        rho = gs.black_hole_mean_density(M)
        mu = (rho / rho_ref) ** (1.0 / 3.0)
        mags.append(mu)
        clean = label.replace("\n", " ").replace("$", "").replace("\\odot", "sun") \
                     .replace("\\", "")
        lines.append(f"  {clean:30s} {rho:.2e}        {mu:.2e}")
    lines += [
        "",
        "READING: an OGU is a horizon-mass black hole, so its MEAN density equals the",
        "cosmic/void density (the R_s = R_H identity that defines a universe). Its",
        "gravity-scaling magnification vs the void is therefore ~1: gravity-scaled",
        "coordinates EQUAL true coordinates at the OGU/void level, and the dilute",
        "scatter (OGUs tiny, round, ~e^{I/4} horizons apart) survives the scaling --",
        "the OGU-level wallpaper is honestly near-empty in EITHER coordinate system.",
        "One level down, astrophysical black holes are ~1e9-1e14x denser than their",
        "host universe, so gravity-scaling blows them up enormously -- THIS is where the",
        "conformal map earns its name, making nested compact structure visible. So the",
        "honest gravity-scaled wallpaper is two-tier: a near-empty dilute scatter of",
        "OGUs, each a richly magnified nest inside. The two-tier behaviour is itself a",
        "consequence of 'a universe is its own Schwarzschild radius'.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "gravity_scale.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: magnification bar per object type.
    labels = [o[0] for o in objs]
    cols = [o[3] for o in objs]
    axL.bar(range(len(objs)), mags, color=cols, alpha=0.85, width=0.6)
    axL.axhline(1.0, color="C3", ls="--", lw=1.4)
    axL.text(len(objs) - 0.5, 1.5, "magnification $=1$\n(gravity-scaled $=$ true)",
             fontsize=8.5, color="C3", ha="right")
    for i, mu in enumerate(mags):
        axL.text(i, mu * 1.5, f"{mu:.0e}" if mu > 10 else f"{mu:.2f}",
                 ha="center", fontsize=8.5)
    axL.set_yscale("log")
    axL.set_xticks(range(len(objs)))
    axL.set_xticklabels(labels, fontsize=8)
    axL.set_ylabel(r"gravity-scaling magnification  $\mu=(\rho/\rho_{\rm ref})^{1/3}$")
    axL.set_title("An OGU is NOT magnified vs its void ($\\mu\\!\\sim\\!1$)\n"
                  "astrophysical black holes are ($\\mu\\!\\sim\\!10^{9\\text{--}14}$)",
                  fontsize=11)
    axL.set_ylim(0.3, 1e16)
    axL.grid(True, axis="y", which="both", alpha=0.2)

    # Panel R: black-hole mean density vs mass, crossing the void/critical line.
    M = np.logspace(30, 54, 300)
    rho = gs.black_hole_mean_density(M)
    axR.plot(M, rho, "C0", lw=2.2, label=r"BH mean density $\propto M^{-2}$")
    axR.axhline(gs.RHO_CRIT, color="C3", ls="--", lw=1.3,
                label=r"cosmic / void density")
    axR.plot([9.2e52], [gs.black_hole_mean_density(9.2e52)], "ko", ms=8)
    axR.annotate("OGU (horizon mass):\nmean density $=$ cosmic\n$\\Rightarrow\\mu=1$",
                 xy=(9.2e52, gs.black_hole_mean_density(9.2e52)),
                 xytext=(3e44, 1e-18), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axR.text(3e31, 1e8, "astrophysical holes:\ndense $\\Rightarrow$ magnified",
             fontsize=8.5, color="0.4")
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"black-hole mass  $M$  [kg]")
    axR.set_ylabel(r"mean density over $R_s$  [kg/m$^3$]")
    axR.set_title("Why: a horizon-mass hole has the cosmic mean density\n"
                  "(the $R_s=R_H$ identity), so an OGU scales to $\\mu=1$",
                  fontsize=11)
    axR.legend(fontsize=8.5, loc="upper right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "gravity_scale.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "gravity_scale.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'gravity_scale.pdf')}")


if __name__ == "__main__":
    main()
