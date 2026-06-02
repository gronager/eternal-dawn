#!/usr/bin/env python3
r"""Tier 2 first result: the primordial spectrum from the Einstein--Cartan bounce.

Left: the tensor power spectrum P_T(k) from modes evolved through the bounce --
flat (scale-invariant) for matter-dominated contraction (w=0), steeply blue for
radiation. Right: the tilt n_T vs the contraction equation of state w, numeric
vs analytic n_T = 3 - |2p-1|, p = 2/(1+3w); matter contraction gives scale
invariance with no inflaton, and the observed slight red tilt needs w just below 0.
Renders figures/pdf/primordial.pdf, writes sims/output/primordial.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import primordial as pr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def w_for_tilt(target_nT):
    """Solve 3-|2p-1| = target for p>1/2, then w = (2/p - 1)/3."""
    p = (3.0 - target_nT + 1.0) / 2.0          # for p>1/2: 3-(2p-1)=target
    return (2.0 / p - 1.0) / 3.0


def main() -> None:
    w_obs = w_for_tilt(pr.N_S_OBS - 1.0)       # contraction w giving the red tilt

    lines = [
        "Primordial spectrum from the Einstein-Cartan bounce (Tier 2 first result)",
        "=" * 62,
        "  modes evolved through the bounce; tensor tilt n_T = dlnP_T/dlnk:",
        "    w (contraction)   p=2/(1+3w)   n_T numeric   n_T analytic",
    ]
    for w in (-0.05, 0.0, 0.05, 1.0 / 3.0):
        p = pr.p_from_w(w)
        lines.append(f"    {w:+.3f}            {p:.3f}        "
                     f"{pr.tilt_numeric(p):+.3f}        {pr.analytic_tilt(p):+.3f}")
    lines += [
        "",
        f"  matter contraction (w=0) -> n_T = 0 EXACTLY: scale-invariant, no inflaton.",
        f"  observed n_s = {pr.N_S_OBS} (slightly red) <-> contraction w ~ {w_obs:+.4f}",
        f"  (just below matter: a mild quintessence-like softening).",
        "",
        "READING: the bounce reproduces the observed near-scale-invariance WITHOUT",
        "inflation -- modes exit the comoving horizon during a matter-dominated",
        "contraction, freeze, and cross the bounce, giving n_T = 3-|2p-1| = 0 at",
        "w=0. What the Einstein-Cartan bounce adds is the non-singular turn itself",
        "(torsion at rho_C), which the matter-bounce scenario otherwise needs exotic",
        "ghost matter to achieve. The precise red tilt n_s=0.965 (a slightly soft,",
        "w<0 contraction) and the tensor-to-scalar ratio r are the remaining Tier-2",
        "work; here the MECHANISM and the scale-invariance are validated.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "primordial.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: P_T(k) for three contractions, normalized at a pivot.
    ks = np.logspace(-2.0, -0.7, 12)
    kpiv = ks[len(ks) // 2]
    cases = [(0.0, "matter $w=0$ (scale-invariant)", "C0", "-"),
             (w_obs, fr"$w={w_obs:+.3f}$ (matches $n_s$)", "C2", "--"),
             (1.0 / 3.0, "radiation $w=1/3$ (blue)", "C3", "-.")]
    for w, lbl, col, ls in cases:
        p = pr.p_from_w(w)
        P = pr.power_spectrum(ks, p)
        Ppiv = pr.power_spectrum([kpiv], p)[0]
        axL.plot(ks, P / Ppiv, col, ls=ls, lw=2.0, label=lbl)
    axL.axhline(1.0, color="0.6", ls=":", lw=0.8)
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"comoving wavenumber  $k$  (bounce units)")
    axL.set_ylabel(r"tensor power  $P_T(k)$  (normalized at pivot)")
    axL.set_title("The bounce sets the primordial spectrum\n"
                  "matter contraction $\\Rightarrow$ scale-invariant, no inflaton",
                  fontsize=11)
    axL.legend(fontsize=8.5, loc="upper left")
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: tilt vs contraction equation of state, near w=0.
    ws = np.linspace(-0.06, 0.10, 200)
    nT_an = np.array([pr.analytic_tilt(pr.p_from_w(w)) for w in ws])
    axR.plot(ws, nT_an, "C0", lw=2.0, label=r"analytic $n_T=3-|2p-1|$")
    ws_pts = np.array([-0.05, -0.02, 0.0, 0.03, 0.06, 0.1])
    axR.plot(ws_pts, [pr.tilt_numeric(pr.p_from_w(w)) for w in ws_pts], "ko",
             ms=5, label="numeric (mode evolution)")
    axR.axhline(0.0, color="0.6", ls=":", lw=0.8)
    axR.axvline(0.0, color="C0", ls="--", lw=1.0, alpha=0.6)
    axR.text(0.002, 1.4, "matter\n(scale-invariant)", fontsize=8, color="C0")
    axR.axhline(pr.N_S_OBS - 1.0, color="C2", ls="--", lw=1.2)
    axR.plot([w_obs], [pr.N_S_OBS - 1.0], "C2*", ms=13)
    axR.annotate(fr"observed $n_s-1={pr.N_S_OBS-1:.3f}$" "\n"
                 fr"$\Rightarrow w\approx{w_obs:+.3f}$",
                 xy=(w_obs, pr.N_S_OBS - 1.0), xytext=(-0.055, -0.55), fontsize=8,
                 color="C2", arrowprops=dict(arrowstyle="->", color="C2"))
    axR.set_xlabel(r"contraction equation of state  $w$")
    axR.set_ylabel(r"tensor tilt  $n_T$")
    axR.set_title("Tilt tracks the contraction\n"
                  "scale-invariance at $w=0$; slight red needs $w\\lesssim0$",
                  fontsize=11)
    axR.set_ylim(-0.8, 2.2)
    axR.legend(fontsize=8.5, loc="upper right")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "primordial.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
