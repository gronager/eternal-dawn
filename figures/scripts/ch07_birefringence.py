#!/usr/bin/env python3
r"""Cosmic birefringence: the parity-odd CMB fingerprint of Eternal Dawn (Ch.7).

Left -- the observable the framework predicts given a rotation beta: the EB and TB
spectra, computed from the REAL CAMB EE/TE transfer via
C_l^{EB}=1/2 sin(4 beta)(EE-BB) and C_l^{TB}=sin(2 beta)TE, for beta=0.1/0.3/0.5 deg.
LambdaCDM predicts these are exactly zero (dashed), so a nonzero EB IS the
birefringence.

Middle -- the published isotropic-beta measurements (Planck, WMAP) against the
LambdaCDM beta=0 line: ~0.3 deg, each ~2.5-3.6 sigma. (They share Planck data, so the
true joint significance is ~3 sigma, not the naive inverse-variance combination.)

Right -- the distinctive, uniquely-ED component: an ANISOTROPIC beta(n) locked to the
inherited parent-spin axis, predicted to coincide with the CMB 'axis of evil' /
galaxy-spin axis (Ch.6). The amplitude is shown as a shape, not a derived number; the
falsifiable content is the AXIS, not the value.

Renders figures/pdf/birefringence.pdf, writes sims/output/birefringence.txt.
Requires CAMB.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import camb_cmb as cc
from cartasis_sims import birefringence as bi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

AXIS_LONLAT = (260.0, 60.0)   # CMB 'axis of evil' direction, Galactic (l, b)


def main():
    if not cc.HAVE_CAMB:
        msg = "CAMB not installed -- run `uv pip install camb` (sims/). Skipping."
        print(msg)
        with open(os.path.join(OUT_DIR, "birefringence.txt"), "w") as f:
            f.write(msg + "\n")
        return

    ell, TT, EE, BB, TE = cc.cmb_polarization_spectra(lmax=2500)
    beta_c, sig_c, signif = bi.combined_beta()

    lines = [
        "Cosmic birefringence: the parity-odd CMB fingerprint of Eternal Dawn",
        "=" * 68,
        "  LambdaCDM: EB = TB = 0 exactly (parity conserved).",
        "  ED: a relic axial-torsion background (Einstein-Cartan) rotates",
        "      polarization by beta, lighting up EB and TB.",
        "",
        "  predicted parity-odd amplitudes from the real CAMB EE/TE:",
        "    beta [deg]   peak |EB|   peak |TB|   [microK^2]",
    ]
    for b in (0.1, 0.3, 0.5):
        eb = bi.eb_peak_amplitude(EE, BB, b)
        tb = float(np.nanmax(np.abs(bi.tb_spectrum(TE, b))))
        lines.append(f"    {b:>6.1f}      {eb:8.3f}    {tb:8.3f}")
    lines += [
        "",
        "  measured isotropic beta (LambdaCDM predicts 0):",
    ]
    for label, b, s, key in bi.BETA_MEASUREMENTS:
        lines.append(f"    {b:.2f} +/- {s:.2f} deg   ({b/s:.1f} sigma)   {label}")
    lines += [
        f"    indicative band: ~{beta_c:.2f} deg "
        "(joint ~3 sigma; the analyses share Planck data)",
        "",
        "READING: ED does not yet predict the VALUE of beta (set by the uncomputed",
        "relic axial-torsion amplitude); it predicts beta != 0 where LambdaCDM forbids",
        "it, and -- distinctively -- an ANISOTROPIC beta(n) whose axis coincides with",
        "the axis-of-evil / galaxy-spin axis (Ch.6). The current ~0.3 deg is",
        "consistent; the decisive ED test is that axis alignment, not the number.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "birefringence.txt"), "w") as f:
        f.write(text + "\n")

    fig = plt.figure(figsize=(15.0, 4.7))
    axA = fig.add_subplot(1, 3, 1)
    axB = fig.add_subplot(1, 3, 2)
    axC = fig.add_subplot(1, 3, 3)

    # ---- A: the EB / TB fingerprint -----------------------------------
    m = (ell >= 2) & (ell <= 2000)
    for b, alpha in [(0.1, 0.4), (0.3, 1.0), (0.5, 0.7)]:
        eb = bi.eb_spectrum(EE, BB, b)
        lw = 2.0 if b == 0.3 else 1.2
        axA.plot(ell[m], eb[m], color="C0", alpha=alpha, lw=lw,
                 label=rf"EB, $\beta={b}^\circ$")
    tb = bi.tb_spectrum(TE, 0.3)
    axA.plot(ell[m], tb[m], color="C1", lw=2.0, label=r"TB, $\beta=0.3^\circ$")
    axA.axhline(0.0, color="k", ls="--", lw=1.0, label=r"$\Lambda$CDM ($=0$)")
    axA.set_xlabel(r"multipole  $\ell$")
    axA.set_ylabel(r"$\mathcal{D}_\ell$  [$\mu$K$^2$]")
    axA.set_title("The parity-odd fingerprint\n(from real CAMB EE/TE)", fontsize=10.5)
    axA.legend(fontsize=8, loc="upper left", ncol=2)
    axA.grid(True, alpha=0.2)

    # ---- B: the measurement forest ------------------------------------
    ys = np.arange(len(bi.BETA_MEASUREMENTS))[::-1]
    for y, (label, b, s, key) in zip(ys, bi.BETA_MEASUREMENTS):
        axB.errorbar(b, y, xerr=s, fmt="C3o", ms=6, capsize=3, lw=1.5)
        axB.text(b, y + 0.18, label.split("(")[-1].rstrip(")"),
                 ha="center", fontsize=7.5, color="0.3")
    axB.axvline(0.0, color="k", ls="--", lw=1.2)
    axB.text(0.005, len(ys) - 0.5, r"$\Lambda$CDM", fontsize=9, color="k")
    axB.axvspan(beta_c - sig_c, beta_c + sig_c, color="C3", alpha=0.12)
    axB.set_xlim(-0.1, 0.6)
    axB.set_ylim(-0.6, len(ys) - 0.2)
    axB.set_yticks([])
    axB.set_xlabel(r"isotropic birefringence  $\beta$  [deg]")
    axB.set_title("Measured: $\\beta\\sim0.3^\\circ$, "
                  "$\\sim3\\sigma$\n($\\Lambda$CDM forbids it)", fontsize=10.5)
    axB.grid(True, axis="x", alpha=0.2)

    # ---- C: the axis-aligned anisotropy (the distinctive part) --------
    lon, lat, bmap = bi.axis_anisotropy_map(0.30, 0.25, axis_lonlat=AXIS_LONLAT)
    im = axC.imshow(bmap, origin="lower", aspect="auto",
                    extent=[lon[0], lon[-1], lat[0], lat[-1]],
                    cmap="RdBu_r")
    axC.plot(AXIS_LONLAT[0] - 360 if AXIS_LONLAT[0] > 180 else AXIS_LONLAT[0],
             AXIS_LONLAT[1], "k*", ms=14)
    axC.annotate("parent-spin axis\n(= axis of evil,\ngalaxy spins)",
                 (AXIS_LONLAT[0] - 360, AXIS_LONLAT[1]),
                 textcoords="offset points", xytext=(8, -28), fontsize=7.5)
    axC.set_xlabel(r"Galactic $\ell$ [deg]")
    axC.set_ylabel(r"Galactic $b$ [deg]")
    axC.set_title("Distinctive: anisotropic $\\beta(\\hat n)$\n"
                  "locked to the inherited axis (schematic)", fontsize=10.5)
    cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.04)
    cb.set_label(r"$\beta(\hat n)$ [deg]", fontsize=8)

    fig.suptitle("Cosmic birefringence: a parity-violating signal "
                 r"$\Lambda$CDM forbids, and the axis that would make it ours",
                 fontsize=12.5, y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "birefringence.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "birefringence.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'birefringence.pdf')}")


if __name__ == "__main__":
    main()
