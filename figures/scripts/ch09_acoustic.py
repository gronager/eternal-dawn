#!/usr/bin/env python3
r"""Toward the acoustic peaks: SCT inherits them from standard recombination.

The bounce gives a near-scale-invariant primordial spectrum (primordial.py); the
acoustic peaks are recombination-era sound waves, shared with LambdaCDM. Left: a
schematic CMB TT spectrum from a scale-invariant primordial tilt x a damped
acoustic transfer, with the observed Planck peak positions marked. Right: the
SCT-computed peak multipoles (from the sound horizon) vs Planck -- they line up,
because it is the same recombination physics on the post-bounce plasma.
Renders figures/pdf/acoustic.pdf, writes sims/output/acoustic.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import acoustic as ac

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

PLANCK_PEAKS = [220.6, 537.5, 810.8, 1120.9, 1444.0]   # observed TT peak l's


def main() -> None:
    theta, lA, rs, dc = ac.acoustic_scale()
    sct_peaks = ac.peak_multipoles(n_peaks=5)

    lines = [
        "Toward the acoustic peaks: SCT inherits them (shared recombination)",
        "=" * 60,
        f"  sound horizon r_s        = {rs:.1f} Mpc   (Planck ~ 144.4)",
        f"  distance to LSS  D_C     = {dc:.0f} Mpc   (Planck ~ 13870)",
        f"  acoustic scale theta_s   = {theta*180/np.pi*60:.2f} arcmin  (Planck ~ 35.7)",
        f"  acoustic multipole l_A   = {lA:.0f}",
        "",
        "  peak multipoles  l_n:   SCT (computed)   vs   Planck (observed)",
    ]
    for i, (s, p) in enumerate(zip(sct_peaks, PLANCK_PEAKS), 1):
        lines.append(f"    peak {i}:            {s:6.0f}             {p:6.0f}")
    lines += [
        "",
        "READING: the bounce supplies the (near scale-invariant) primordial spectrum;",
        "the acoustic peaks are pre-recombination baryon-photon sound waves, which",
        "SCT shares entirely with LambdaCDM -- same hot Big Bang after the bounce,",
        "same recombination at z~1090, same sound horizon. So the peak SERIES is",
        "inherited, and the computed positions match Planck. The 'CMB = parent",
        "Hawking radiation' identification concerns the ORIGIN of the thermal plasma",
        "at the bounce surface, not a replacement for recombination; the peaks come",
        "out either way. (The schematic C_l is illustrative; a Boltzmann computation",
        "is the next step, but the acoustic SCALE -- the hard, decisive number -- is",
        "reproduced.)",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "acoustic.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: schematic TT spectrum with Planck peak positions.
    ell = np.linspace(2, 2200, 1500)
    Dl = ac.schematic_cl(ell)
    Dl = Dl / Dl.max()
    axL.plot(ell, Dl, "C0", lw=2.0, label="schematic $D_\\ell$ (bounce primordial)")
    for i, lp in enumerate(PLANCK_PEAKS):
        axL.axvline(lp, color="C3", ls=":", lw=1.1,
                    label="Planck peak $\\ell$" if i == 0 else None)
    axL.set_xlabel(r"multipole  $\ell$")
    axL.set_ylabel(r"$\mathcal{D}_\ell=\ell(\ell+1)C_\ell/2\pi$  (normalized)")
    axL.set_title("Acoustic peaks ride on the bounce's scale-invariant input\n"
                  "(schematic transfer; peak positions are the physics)", fontsize=10.5)
    axL.legend(fontsize=8.5, loc="upper right")
    axL.grid(True, alpha=0.2)

    # Panel R: SCT-computed peaks vs Planck (1:1).
    axR.plot(PLANCK_PEAKS, sct_peaks, "C0o", ms=8)
    lim = [0, 1600]
    axR.plot(lim, lim, "0.5", ls="--", lw=1.0, label="1:1")
    for s, p, i in zip(sct_peaks, PLANCK_PEAKS, range(1, 6)):
        axR.annotate(f"peak {i}", (p, s), xytext=(6, -10),
                     textcoords="offset points", fontsize=8)
    axR.set_xlim(*lim)
    axR.set_ylim(*lim)
    axR.set_xlabel(r"Planck observed peak  $\ell_n$")
    axR.set_ylabel(r"SCT computed peak  $\ell_n=(n-\frac{1}{4})\pi/\theta_s$")
    axR.set_title(f"Peak positions match\n"
                  f"$r_s={rs:.0f}$ Mpc, $\\theta_s={theta*180/np.pi*60:.1f}'$"
                  f" (shared recombination)", fontsize=10.5)
    axR.legend(fontsize=8.5, loc="upper left")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "acoustic.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
