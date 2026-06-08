#!/usr/bin/env python3
r"""The torsiton tower as a finite Woods-Saxon bag: three rungs, the cap, and the mass overlap.

Left: the Woods-Saxon well V(r) tuned to bind EXACTLY three s-wave rungs (the generations) -- the
ground (0 nodes), and two radial excitations (1 and 2 nodes) -- each wavefunction drawn riding on its
energy level. The condensate (the bag-filling chiral substrate) is the shaded profile. A fourth state
would sit above the dissociation threshold (V=0): the cap. Right: the configurational mass of each
rung, its overlap with the condensate (m_n = int u_n^2 sigma dr) -- the mass mechanism, and how it
falls as the rung spreads and develops nodes. Renders figures/pdf/woods_saxon_torsiton.pdf.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import woods_saxon as ws

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    V0, R0, a = 5.0, 3.0, 0.5
    r, states = ws.ws_spectrum(V0, R0, a, rmax=14.0, N=1400)
    V = ws.ws_potential(r, V0, R0, a)
    cond = 1.0 / (1.0 + np.exp((r - R0) / a))                 # condensate: the bag-filling substrate
    masses = ws.overlap_masses(r, states, cond)
    masses = masses / masses.max()                            # normalise (heaviest rung = 1)
    nodes = [ws.interior_nodes(u) for _, u in states]
    labels = ["gen I", "gen II", "gen III"]
    cols = ["C3", "C0", "C2"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 4.6), gridspec_kw={"width_ratios": [2.0, 1.0]})

    # --- left: the well, the rungs riding on their levels, the condensate ---
    axL.fill_between(r, V, 0.2, color="0.92", zorder=0)        # the well body
    axL.plot(r, V, color="0.35", lw=1.6)
    axL.fill_between(r, -V0 - 0.3, -V0 - 0.3 + 1.4 * cond, color="C1", alpha=0.18, zorder=0)
    axL.plot(r, -V0 - 0.3 + 1.4 * cond, color="C1", lw=1.4, label="condensate $\\sigma(r)$")
    axL.axhline(0.0, color="0.5", ls=":", lw=1)               # dissociation threshold
    axL.text(11.6, 0.08, "threshold (4th rung\nwould dissociate)", fontsize=8, color="0.4", va="bottom")
    sc = 1.5
    for n, (E, u) in enumerate(states):
        axL.axhline(E, xmin=0.02, xmax=0.62, color=cols[n], ls="--", lw=0.8, alpha=0.6)
        axL.plot(r, E + sc * u, color=cols[n], lw=2.0)
        axL.fill_between(r, E, E + sc * u, color=cols[n], alpha=0.12)
        axL.text(8.6, E, f"{labels[n]}  ({nodes[n]} node{'s' if nodes[n]!=1 else ''})",
                 color=cols[n], fontsize=9, va="center")
    axL.set_xlim(0, 13); axL.set_ylim(-V0 - 0.6, 1.4)
    axL.set_xlabel("$r$"); axL.set_ylabel("energy  /  $V(r)$")
    axL.set_title("the torsiton bag: three rungs, capped", fontsize=11)
    axL.legend(loc="lower right", fontsize=8.5)

    # --- right: the overlap mass per generation (the mass mechanism) ---
    x = np.arange(len(states))
    axR.bar(x, masses, color=cols[:len(states)], alpha=0.85)
    for n in range(len(states)):
        axR.text(n, masses[n] + 0.02, f"{masses[n]:.2f}", ha="center", fontsize=9)
    axR.set_xticks(x); axR.set_xticklabels(labels[:len(states)])
    axR.set_ylabel(r"overlap mass  $m_n=\int u_n^2\,\sigma\,dr$  (norm.)")
    axR.set_title("mass = overlap with $\\sigma$", fontsize=11)
    axR.set_ylim(0, 1.15)
    axR.text(0.5, 0.5, "the rung spreads &\ngains nodes $\\Rightarrow$\nless overlap $\\Rightarrow$\nlighter",
             transform=axR.transAxes, fontsize=8.5, color="0.35", ha="center", va="center")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "woods_saxon_torsiton.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  V0={V0} R0={R0} a={a}: {len(states)} bound rungs  E={np.round([e for e,_ in states],3)}")
    print(f"  interior nodes: {nodes}   overlap masses (norm): {np.round(masses,3)}")


if __name__ == "__main__":
    main()
