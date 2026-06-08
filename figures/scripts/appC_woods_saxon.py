#!/usr/bin/env python3
r"""The torsiton tower as a finite Woods-Saxon bag: three rungs, the cap, and the mass overlap.

Left: the Woods-Saxon well V(r) tuned to bind EXACTLY three s-wave rungs (the generations) -- the
ground (0 nodes), and two radial excitations (1 and 2 nodes) -- each wavefunction drawn riding on its
energy level. Two complementary profiles fill the picture: the condensate sigma(r) (large INSIDE the
bag, the chiral-restored floor) and the local dynamical mass M(r) = M_vac*(1 - sigma/sigma_max) (small
in the core, large in the VACUUM outside). A fourth state would sit above the dissociation threshold
(V=0): the cap. Right: the configurational mass of each rung is its overlap with M(r), NOT with sigma
-- Eq. configmass weights the scalar density by the LOCAL MASS, which is large in the vacuum. So the
mass RISES with generation: the ground rung is tucked in the core where M is small (the electron, the
lightest), while higher rungs spread out and reach into the vacuum where M is large (heavier). This is
the physical ordering. Renders figures/pdf/woods_saxon_torsiton.pdf.
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
    V0, R0, a = 1.0, 5.0, 1.0                                 # tuned: 3 rungs AND ground rms ~ torsiton
    r, states = ws.ws_spectrum(V0, R0, a, rmax=16.0, N=1600)
    V = ws.ws_potential(r, V0, R0, a)
    cond = 1.0 / (1.0 + np.exp((r - R0) / a))                 # condensate sigma(r): large INSIDE the bag
    mass_r = 1.0 - cond                                       # local dynamical mass M(r): large in VACUUM
    masses = ws.overlap_masses(r, states, mass_r)            # Eq. configmass: weight the density by M(r)
    masses = masses / masses[0]                               # normalise to the ground rung (electron=1)
    nodes = [ws.interior_nodes(u) for _, u in states]
    labels = ["gen I", "gen II", "gen III"]
    cols = ["C3", "C0", "C2"]

    # the ACTUAL torsiton ground state (self-consistent well, lattice-confirmed to bind), for overlay
    from cartasis_sims import self_consistent as scc
    from cartasis_sims.self_consistent import _solve_levels
    o = scc.solve_soliton(m0=1.0, g=4.0, m_sigma=1.5, n_fermions=1, n_levels=4, R=16.0, N=1600,
                          iters=900, mix=0.12)
    rt = o["r"]
    _, ut = _solve_levels(o["M"], rt, rt[1] - rt[0], 4)
    u_t = ut[:, 0]
    u_t = u_t / np.max(np.abs(u_t))                           # normalise to unit peak for shape overlay

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 4.6), gridspec_kw={"width_ratios": [2.0, 1.0]})

    # --- left: the well, the rungs riding on their levels, the condensate ---
    axL.fill_between(r, V, 0.15, color="0.92", zorder=0)       # the well body
    axL.plot(r, V, color="0.35", lw=1.6)
    base = -V0 - 0.5
    axL.fill_between(r, base, base + 0.45 * cond, color="C1", alpha=0.18, zorder=0)
    axL.plot(r, base + 0.45 * cond, color="C1", lw=1.4, label="condensate $\\sigma(r)$ (in the bag)")
    axL.plot(r, base + 0.45 * mass_r, color="C4", lw=1.4, ls="-.",
             label="local mass $M(r)$ (in the vacuum)")
    axL.axhline(0.0, color="0.5", ls=":", lw=1)               # dissociation threshold
    axL.text(11.4, 0.02, "threshold\n(4th dissociates)", fontsize=8, color="0.4", va="bottom")
    sc = 0.30
    E0 = states[0][0]
    u0 = states[0][1] / np.max(np.abs(states[0][1]))          # WS ground, unit peak
    for n, (E, u) in enumerate(states):
        un = u / np.max(np.abs(states[0][1]))                 # common scale (WS ground peak = 1)
        axL.axhline(E, xmin=0.02, xmax=0.60, color=cols[n], ls="--", lw=0.7, alpha=0.5)
        axL.plot(r, E + sc * un, color=cols[n], lw=2.0)
        axL.fill_between(r, E, E + sc * un, color=cols[n], alpha=0.10)
        axL.text(9.4, E, f"{labels[n]}  ({nodes[n]} node{'s' if nodes[n]!=1 else ''})",
                 color=cols[n], fontsize=9, va="center")
    # overlay the ACTUAL torsiton ground state on the gen-I level (sign-matched)
    if np.sum(u0 * np.interp(r, rt, u_t)) < 0:
        u_t = -u_t
    axL.plot(r, E0 + sc * np.interp(r, rt, u_t), color="k", lw=1.7, ls=(0, (4, 2)),
             label="actual torsiton $\\psi_0$ (self-consistent)")
    axL.set_xlim(0, 13); axL.set_ylim(-V0 - 0.7, 0.5)
    axL.set_xlabel("$r$"); axL.set_ylabel("energy  /  $V(r)$")
    axL.set_title("the torsiton bag: three rungs, capped (shape-matched to $\\psi_0$)", fontsize=10.5)
    axL.legend(loc="lower right", fontsize=8)

    # --- right: the overlap mass per generation (the mass mechanism) ---
    x = np.arange(len(states))
    axR.bar(x, masses, color=cols[:len(states)], alpha=0.85)
    for n in range(len(states)):
        axR.text(n, masses[n] + 0.04 * masses.max(), f"{masses[n]:.2f}", ha="center", fontsize=9)
    axR.set_xticks(x); axR.set_xticklabels(labels[:len(states)])
    axR.set_ylabel(r"configurational mass  $m_n=\int u_n^2\,M(r)\,dr$  (norm.)")
    axR.set_title("mass = overlap with $M(r)$", fontsize=11)
    axR.set_ylim(0, 1.18 * masses.max())
    axR.text(0.30, 0.78, "higher rungs spread\ninto the vacuum, where\n$M(r)$ is large $\\Rightarrow$\nheavier",
             transform=axR.transAxes, fontsize=8.0, color="0.35", ha="center", va="center")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "woods_saxon_torsiton.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  V0={V0} R0={R0} a={a}: {len(states)} bound rungs  E={np.round([e for e,_ in states],3)}")
    print(f"  interior nodes: {nodes}   overlap masses (norm): {np.round(masses,3)}")


if __name__ == "__main__":
    main()
