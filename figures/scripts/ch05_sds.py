#!/usr/bin/env python3
r"""Why an OGU evaporates though it looks "in equilibrium": the two SdS horizons.

Left: Schwarzschild-de Sitter has a black-hole horizon r_b and a cosmological horizon
r_c; they merge at the Nariai mass. Right: the two temperatures -- T_b is ALWAYS above
T_c (equal only at Nariai), so net radiation always flows outward and the OGU
evaporates. The "equilibrium" was a near-Nariai lukewarm plateau (T_b - T_c -> 0),
not true balance, and the mass cap is geometric (no bigger BH fits), not thermal.
Renders figures/pdf/sds.pdf, writes sims/output/sds.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import sds

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    fracs = np.linspace(0.002, 0.99995, 500)
    mus = fracs * sds.MU_NARIAI
    rb = np.array([sds.sds_horizons(m)[0] for m in mus])
    rc = np.array([sds.sds_horizons(m)[1] for m in mus])
    Tb = np.array([sds.sds_temperatures(m)[0] for m in mus])
    Tc = np.array([sds.sds_temperatures(m)[1] for m in mus])

    lines = [
        "Schwarzschild-de Sitter: why the OGU evaporates (no true equilibrium)",
        "=" * 62,
        f"  Nariai mass mu_N = 1/(3 sqrt3) = {sds.MU_NARIAI:.4f} c^3/(G H_Lambda)",
        f"  at Nariai the two horizons merge: r_b = r_c = 1/sqrt3 = {sds.R_NARIAI:.4f}",
        "",
        "    M/M_Nariai   T_b/T_c   net outflow?",
    ]
    for f in (0.2, 0.5, 0.95, 0.999):
        mu = f * sds.MU_NARIAI
        lines.append(f"    {f:.3f}        {sds.temperature_ratio(mu):.3f}     "
                     f"{sds.net_outflow(mu)}")
    lines += [
        "",
        "RESOLUTION: T_b > T_c for EVERY sub-Nariai mass (equal only at the unstable",
        "Nariai point). The black hole is always hotter than its own cosmological",
        "horizon, so net radiation always flows outward and the OGU always evaporates.",
        "The earlier 'reaches equilibrium, growth stops' was wrong twice over:",
        "  * the mass cap is GEOMETRIC (no larger black hole fits in this de Sitter),",
        "    not a thermal balance;",
        "  * near Nariai T_b - T_c -> 0, so evaporation becomes very slow -- a lukewarm",
        "    near-degenerate plateau that LOOKS like equilibrium but is not.",
        "So the OGU caps at ~Nariai (geometry) and then drains over ~1e142 s (residual",
        "T_b > T_c). No contradiction.",
        "",
        "DOES TIME INSIDE STOP WHEN THE MEMBRANE EVAPORATES? No. The membrane is the",
        "baby's PAST horizon (its Big-Bang surface), not its future; the baby's own",
        "time runs to its own de Sitter infinity, eternally. Exterior evaporation and",
        "interior eternity are DUALS across the null horizon: an infinite interior",
        "proper time maps to the finite exterior horizon lifetime. By the time the",
        "parent horizon has evaporated, the baby has already (on its own clock) reached",
        "its de Sitter future = become a flat heat-dead void hosting its own sub-OGUs.",
        "The parent's evaporation is just the exterior face of the interior reaching",
        "conformal infinity: the umbilical channel pinches off, but no clock stops and",
        "no nested structure is destroyed -- the inner universes were always inside the",
        "baby's own far-future void, which simply becomes a free-standing piece of the",
        "eternal flat substrate. 'No discontinuities' all the way down.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "sds.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: the two horizon radii vs mass, merging at Nariai.
    axL.plot(fracs, rb, "C0", lw=2.2, label=r"black-hole horizon $r_b$")
    axL.plot(fracs, rc, "C3", lw=2.2, label=r"cosmological horizon $r_c$")
    axL.fill_between(fracs, rb, rc, color="0.5", alpha=0.12)
    axL.text(0.3, 0.78, "the universe\n(between the horizons)", fontsize=8.5,
             color="0.35")
    axL.plot([1.0], [sds.R_NARIAI], "ko", ms=7)
    axL.annotate("Nariai: horizons merge\n(geometric max mass)",
                 xy=(1.0, sds.R_NARIAI), xytext=(0.45, 0.30), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.set_xlabel(r"OGU mass  $M/M_{\rm Nariai}$")
    axL.set_ylabel(r"horizon radius  [$c/H_\Lambda$]")
    axL.set_title("Schwarzschild--de Sitter has two horizons\n"
                  "they merge at the Nariai (maximum) mass", fontsize=11)
    axL.legend(fontsize=8.5, loc="center left")
    axL.grid(True, alpha=0.2)

    # Panel R: the two temperatures -- T_b always above T_c.
    axR.plot(fracs, Tb, "C0", lw=2.2, label=r"$T_b$ (black-hole horizon)")
    axR.plot(fracs, Tc, "C3", lw=2.2, label=r"$T_c$ (cosmological horizon)")
    axR.fill_between(fracs, Tc, Tb, color="C0", alpha=0.12)
    axR.text(0.35, 0.45, "$T_b>T_c$ always\n$\\Rightarrow$ net outflow\n"
             "$\\Rightarrow$ evaporates", fontsize=9, color="C0")
    axR.plot([1.0], [Tb[-1]], "ko", ms=6)
    axR.annotate("Nariai: $T_b=T_c$\n(the only 'equilibrium',\nunstable, never reached)",
                 xy=(1.0, Tb[-1]), xytext=(0.5, 0.9), fontsize=8,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axR.set_xlabel(r"OGU mass  $M/M_{\rm Nariai}$")
    axR.set_ylabel(r"horizon temperature  [$H_\Lambda/2\pi$]")
    axR.set_title("The black hole is always hotter than its void\n"
                  "so it drains outward -- slowly near Nariai (lukewarm)", fontsize=11)
    axR.set_yscale("log")
    axR.set_ylim(1e-2, 5.0)
    axR.legend(fontsize=8.5, loc="upper right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "sds.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
