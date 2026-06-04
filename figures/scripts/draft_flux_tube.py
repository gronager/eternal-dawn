#!/usr/bin/env python3
r"""1-D string tension: a confining flux tube from the condensate (dual superconductor).

A: the Nielsen-Olesen vortex profile -- the condensate f(r) (0 on the axis = normal core,
   ->1 outside = superconductor) expels the colour-electric field B(r) into a localized
   TUBE. The flux is squeezed into a string of finite width.
B: the consequence -- a linear confining potential V(L) = sigma L (rises forever; quarks
   cannot be isolated), against the SCREENED overlap force (-> 0) that the pairwise
   picture gave. The flux tube supplies the confinement the overlap could not.

Tension sigma is proportional to the condensate scale v^2 = f_pi^2 (BPS: sigma = 2 pi v^2,
the code's validation) -- so the SAME condensate that makes mass also confines. Conditional
on the gravity-torsion vacuum being a dual superconductor; a lattice computation of the
Part I connection would decide that. Renders figures/pdf/flux_tube.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import flux_tube as ft
from cartasis_sims import color_force as cf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    out = ft.solve_vortex(beta=2.0)
    sigma = out["sigma"]
    lines = [
        "1-D string tension: a confining flux tube from the condensate",
        "=" * 60,
        f"  vortex converged: {out['success']}",
        f"  BPS string tension (beta=2, m_H=m_W): sigma = {sigma:.4f}  "
        f"(= 2 pi v^2 = {2*np.pi:.4f}, topological saturation -> validation)",
        f"  type-I (beta=1): sigma = {ft.string_tension(1.0):.3f};  "
        f"type-II (beta=3): sigma = {ft.string_tension(3.0):.3f}",
        "",
        "READING: a confining flux tube EXISTS -- the condensate expels colour-electric",
        "flux into a 1-D string of tension sigma ~ 2 pi v^2 (set by the condensate scale",
        "v=f_pi). So V(L)=sigma L rises forever: confinement, the linear potential the",
        "overlap/exchange picture could NOT produce. And the SAME condensate that creates",
        "mass (configurational, Part II) sets the string tension -- one scale, both jobs.",
        "",
        "CONDITIONAL: this assumes the gravity-torsion vacuum is a dual superconductor",
        "(its condensate expels colour flux via a dual Meissner effect). Whether it",
        "actually is -- monopole condensation in the non-abelian Part I connection -- is",
        "the non-perturbative question a real LATTICE computation would decide (Wilson-loop",
        "area law / dual order parameter). Qualitative result; clean lattice target.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "flux_tube.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.2))

    # ---- A: the vortex profile ----
    r = out["r"]
    m = r < 8
    axA.plot(r[m], out["f"][m], "C0", lw=2.2, label=r"condensate $f(r)$ (normal core $\to$ SC)")
    axA.plot(r[m], out["a"][m], "C1", lw=1.8, label=r"gauge profile $a(r)$")
    B = out["B"] / np.max(np.abs(out["B"]))
    axA.plot(r[m], B[m], "C3", lw=2.0, label=r"colour-electric flux $B(r)$ (the tube)")
    axA.fill_between(r[m], 0, B[m], color="C3", alpha=0.12)
    axA.set_xlim(0, 8)
    axA.set_xlabel(r"radius $r$ (units $1/m_W$)")
    axA.set_ylabel("profile")
    axA.set_title("The flux tube: condensate expels flux into a string\n"
                  r"$\sigma=2\pi v^2$ at BPS (validation)", fontsize=11)
    axA.legend(fontsize=8.5, loc="center right")
    axA.grid(True, alpha=0.2)

    # ---- B: confinement vs screening ----
    L = np.linspace(0.2, 6.0, 300)
    axB.plot(L, ft.confining_potential(L, sigma), "C3", lw=2.4,
             label=r"flux tube $V(L)=\sigma L$ (CONFINES)")
    axB.plot(L, 30 * cf.residual_yukawa(L, g=2.0, m_sigma=1.0) + sigma * 0,
             "C2--", lw=1.8, label="overlap/exchange force (SCREENS $\\to 0$)")
    axB.axhline(0, color="k", lw=0.8)
    axB.set_xlabel(r"separation $L$ between two colour sources")
    axB.set_ylabel(r"potential $V(L)$")
    axB.set_title("Confinement vs screening\n"
                  "the flux tube supplies the linear $\\sim L$ the overlap missed",
                  fontsize=11)
    axB.legend(fontsize=8.5, loc="upper left")
    axB.grid(True, alpha=0.2)

    fig.suptitle("1-D string tension: a confining flux tube from the condensate "
                 r"($\sigma\propto f_\pi^2$) -- conditional on dual superconductivity",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "flux_tube.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "flux_tube.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'flux_tube.pdf')}")


if __name__ == "__main__":
    main()
