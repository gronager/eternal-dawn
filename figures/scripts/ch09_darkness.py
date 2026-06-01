#!/usr/bin/env python3
r"""Darkness is the harbinger of light: how the dark sector sculpts the CMB.

In SCT the dark sector IS the parent's gravity felt before its light (lookback.py):
dark matter = the parent's accumulated projected mass, dark energy = the time-
derivative of its growth. That same darkness shapes the oldest light we see -- the
CMB -- in two distinct, computable ways (CAMB):

  * DARK MATTER (omega_c) sets the peak HEIGHTS. More cold dark matter deepens the
    potential wells the photon-baryon fluid falls into and damps the radiation
    driving, lowering the first/second peaks and lifting the third relative to them.
    The peak-height envelope is dark matter's fingerprint.
  * DARK ENERGY (Omega_Lambda) sets the peak POSITIONS. It fixes the distance to last
    scattering, hence the angle the sound horizon subtends, hence where every peak
    sits in multipole. More dark energy slides the whole peak series.

So the darkness that fills our metric NOW (felt, not yet seen) is exactly what
sculpted the light emitted at our past horizon THEN: dark matter draws the heights,
dark energy places them. Darkness is the harbinger of light. Renders
figures/pdf/darkness.pdf, writes sims/output/darkness.txt. Requires CAMB.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import camb_cmb as cc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    if not cc.HAVE_CAMB:
        msg = "CAMB not installed -- run `uv pip install camb` (sims/). Skipping."
        print(msg)
        with open(os.path.join(OUT_DIR, "darkness.txt"), "w") as f:
            f.write(msg + "\n")
        return

    ell0, TT0 = cc.cmb_tt_spectrum(lmax=2200)
    base = cc.peak_table()

    lines = [
        "Darkness is the harbinger of light: the dark sector sculpts the CMB",
        "=" * 64,
        "  In SCT the dark sector is the parent's gravity felt before its light",
        "  (lookback.py). That same darkness shapes the oldest light, the CMB:",
        "",
        "  DARK MATTER (omega_c) -> peak HEIGHTS (radiation driving / potential wells):",
        "    omch2   P1      P2      P3     P1/P3",
    ]
    for omc in (0.09, 0.12, 0.15):
        p = cc.peak_table(omch2=omc)
        lines.append(f"    {omc:.2f}   {p[0][1]:5.0f}  {p[1][1]:5.0f}  {p[2][1]:5.0f}"
                     f"   {p[0][1]/p[2][1]:.2f}")
    lines += [
        "",
        "  DARK ENERGY (Omega_Lambda, via H0) -> peak POSITIONS (distance to LSS):",
        "    H0     l1    l2    l3",
    ]
    for H0 in (60.0, 67.36, 75.0):
        p = cc.peak_table(H0=H0)
        lines.append(f"    {H0:5.1f}  {p[0][0]:4d}  {p[1][0]:4d}  {p[2][0]:4d}")
    lines += [
        "",
        "READING: dark matter draws the peak HEIGHTS, dark energy places their",
        "POSITIONS -- between them the dark sector writes both axes of the CMB. In SCT",
        "both are the same thing: the parent's mass-energy reaching us through the",
        "membrane as gravity ahead of light (dark matter ~ its accumulated mass, dark",
        "energy ~ the rate of its growth). So the invisible darkness filling our metric",
        "NOW is exactly what sculpted the visible light emitted at our past horizon",
        "THEN. CMB and dark sector are one phenomenon seen at two times: darkness is",
        "the harbinger of light.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "darkness.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    m = (ell0 >= 2) & (ell0 <= 2200)

    # Panel L: dark matter -> heights. Vary omega_c, hold geometry.
    grays = ["0.0", "0.45", "0.7"]
    for omc, g, ls in zip((0.09, 0.12, 0.15), grays, ("-", "-", "-")):
        ell, TT = cc.cmb_tt_spectrum(lmax=2200, omch2=omc)
        lab = f"$\\omega_c={omc}$" + ("  (ours)" if abs(omc - 0.12) < 1e-9 else "")
        axL.plot(ell[m], TT[m], color=("C0" if abs(omc-0.12)<1e-9 else g),
                 lw=2.2 if abs(omc-0.12)<1e-9 else 1.4, label=lab)
    axL.set_xlabel(r"multipole  $\ell$")
    axL.set_ylabel(r"$\mathcal{D}_\ell$  [$\mu$K$^2$]")
    axL.set_title("Dark matter draws the peak HEIGHTS\n"
                  "(more $\\omega_c$ $\\Rightarrow$ less driving, lower 1st peak)",
                  fontsize=11)
    axL.set_xlim(0, 1500)
    axL.set_ylim(0, 7000)
    axL.legend(fontsize=9, title="cold dark matter", loc="upper right")
    axL.grid(True, alpha=0.2)

    # Panel R: dark energy -> positions. Vary H0 (Omega_Lambda), hold phys densities.
    for H0, g in zip((60.0, 67.36, 75.0), grays):
        ell, TT = cc.cmb_tt_spectrum(lmax=2200, H0=H0)
        OL = 1 - (0.02237 + 0.12) / (H0 / 100) ** 2
        lab = f"$H_0={H0:g}$, $\\Omega_\\Lambda\\approx{OL:.2f}$" + (
            "  (ours)" if abs(H0 - 67.36) < 1e-6 else "")
        axR.plot(ell[m], TT[m], color=("C0" if abs(H0-67.36)<1e-6 else g),
                 lw=2.2 if abs(H0-67.36)<1e-6 else 1.4, label=lab)
    axR.set_xlabel(r"multipole  $\ell$")
    axR.set_ylabel(r"$\mathcal{D}_\ell$  [$\mu$K$^2$]")
    axR.set_title("Dark energy places the peak POSITIONS\n"
                  "(more $\\Omega_\\Lambda$ $\\Rightarrow$ peaks slide to lower $\\ell$)",
                  fontsize=11)
    axR.set_xlim(0, 1500)
    axR.set_ylim(0, 7000)
    axR.legend(fontsize=8.5, title="dark energy", loc="upper right")
    axR.grid(True, alpha=0.2)

    fig.suptitle("Darkness is the harbinger of light: the dark sector sculpts the CMB"
                 "  (dark matter = heights, dark energy = positions)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "darkness.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(PDF_DIR, "darkness.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'darkness.pdf')}")


if __name__ == "__main__":
    main()
