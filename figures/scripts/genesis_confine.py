#!/usr/bin/env python3
r"""Horizon 3a: confinement -- why we see baryons and leptons, never a free quark.

Left: the confining string V(r)=sigma r (the lattice's measured potential) plus the short-range core
that gives a hadron its finite size r*. Middle: the species energetics vs box size -- a single quark's
energy GROWS without bound (the colour string -> confined, no free state), while lepton, meson, and
baryon are finite and box-independent (free). Right: the N-body binding -- coloured quarks (+/-) pulled
into colour-singlet hadrons by the string, while colourless leptons diffuse free. Renders
figures/pdf/genesis_confine.pdf.

HONEST: a U(1)-colour proxy of confinement (the exact-in-2D linear string), the mechanism made
tractable; the full SU(3) flux-tube baryon and the electroweak lepton sector are scoped in HORIZON.md.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_confine as gc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

SIGMA, CORE = 0.08, 4.0


def main():
    fig, (axV, axE, axB) = plt.subplots(1, 3, figsize=(13.4, 4.3))

    # --- left: the confining string + the bound hadron size ---
    r = np.linspace(0.4, 18, 400)
    rstar = gc.meson_size(SIGMA, CORE)
    axV.plot(r, SIGMA * r, "C3-", lw=1.6, label=r"colour string  $V=\sigma r$")
    axV.plot(r, gc.string_potential(r, SIGMA, CORE), "C0-", lw=2.0, label=r"hadron  $\sigma r + c/r$")
    axV.axvline(rstar, color="0.5", ls=":", lw=1.2)
    axV.text(rstar, 0.2, f" $r^*$={rstar:.1f}\n(hadron size)", fontsize=8.5, color="0.4", va="bottom")
    axV.set_xlabel("quark separation $r$"); axV.set_ylabel("energy")
    axV.set_title("the string confines;\nthe core sets a finite size", fontsize=10)
    axV.legend(fontsize=8.5, loc="upper left"); axV.grid(alpha=0.2); axV.set_ylim(0, SIGMA * 18)

    # --- middle: species energetics vs box -- quark confined, others free ---
    boxes = np.linspace(5, 120, 60)
    style = {"lepton": ("C2", "lepton (free)"), "quark": ("C3", "quark (CONFINED)"),
             "meson": ("C0", "meson (free)"), "baryon": ("C1", "baryon (free)")}
    for k, (c, lab) in style.items():
        E = [gc.species_energy(k, SIGMA, box=b, core=CORE) for b in boxes]
        axE.plot(boxes, E, color=c, lw=2.0, label=lab)
    axE.set_xlabel("box size (try to free the colour $\\rightarrow$)")
    axE.set_ylabel("species energy   ($M_0=1$)")
    axE.set_title("a quark has NO free state\n(its energy grows; the rest are flat)", fontsize=10)
    axE.legend(fontsize=8.5, loc="upper left"); axE.grid(alpha=0.2)

    # --- right: N-body binding -- quarks bound, leptons free ---
    res = gc.run_binding(n_quark_pairs=10, n_lepton=8, L=40.0, steps=4000, sigma=SIGMA, core=CORE, seed=3)
    x, q = res["x"], res["q"]
    for charge, col, mk, lab in [(1, "C3", "o", "quark"), (-1, "C0", "o", "antiquark"),
                                 (0, "0.6", "s", "lepton (free)")]:
        m = q == charge
        axB.scatter(x[m, 0], x[m, 1], c=col, marker=mk, s=42 if charge else 30,
                    edgecolors="k", linewidths=0.4, label=lab)
    axB.set_xlim(0, 40); axB.set_ylim(0, 40); axB.set_xticks([]); axB.set_yticks([])
    axB.set_title(f"quarks bound into hadrons ({res['bound_fraction']*100:.0f}%),\nleptons free",
                  fontsize=10)
    axB.legend(fontsize=8, loc="upper right")

    fig.suptitle("Confinement: the strong-torsion string binds quarks into baryons & mesons, "
                 "leaving leptons free  —  why no bare quark", fontsize=11.5, y=1.02)
    fig.tight_layout()
    out = os.path.join(PDF_DIR, "genesis_confine.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  bound fraction {res['bound_fraction']:.2f}, r*={res['r_star']:.1f}, sigma={SIGMA}")


if __name__ == "__main__":
    main()
