#!/usr/bin/env python3
r"""3D genesis: baryons (torsitons) condensing from the cooling SU(2) chiral fluid (Skyrme model).

Top row: max-intensity projections (along z) of the baryon density |B(x)| as the hot, disordered
post-bounce chiral fluid is cooled through the transition -- the baryons crystallise out as bright
localized knots. Bottom: the condensation -- energy and baryon content sum|B| collapse as the fluid
orders and the relic baryon set freezes. This is the real nuclear Skyrme model in 3D, so the knots
are literally baryons (B=1 each); with kappa,mu2 pinned to the lattice scale the knot energy is a
torsiton mass (~30%) and the relic density is the genesis abundance. Renders genesis_quench3d.pdf.

HONEST: classical stochastic EFT (formation + mass-scale + abundance), not the exact quantum spectrum.
Modest grid (numpy); for production swap to cupy/jax on a GPU. The 3D picture, runnable.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_quench3d as g3

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    res = g3.run_quench3d(L=48, steps=2200, T0=0.85, T1=0.0, seed=5, record=24)
    frames = res["frames"]
    idx = [int(0.10 * len(frames)), int(0.35 * len(frames)), int(0.65 * len(frames)), len(frames) - 1]

    fig = plt.figure(figsize=(12.6, 5.0))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.35, 1.0], hspace=0.36, wspace=0.16)
    labels = ["hot tally\n(disordered)", "cooling\n(domains)", "condensing\n(baryons forming)",
              "frozen relic\n(baryons)"]
    vmax = np.percentile(np.abs(frames[idx[2]][1]), 99.9) or 1.0
    for k, i in enumerate(idx):
        s, b = frames[i]
        proj = np.max(np.abs(b), axis=2)                    # max-intensity projection along z
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(proj.T, origin="lower", cmap="inferno", vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{labels[k]}\nT={res['T'][i]:.2f}", fontsize=9.5)
        if k == 0:
            ax.set_ylabel("baryon density $|B(x)|$\n(max-projection)", fontsize=9)

    axc = fig.add_subplot(gs[1, :2])
    t = np.arange(len(res["content"]))
    axc.plot(t, res["content"] / res["content"].max(), "o-", color="C3", ms=3,
             label="baryon content $\\sum|B|$")
    axc.plot(t, res["energy"] / res["energy"].max(), "s-", color="0.5", ms=3, label="energy")
    axc.plot(t, res["T"] / res["T"].max(), ":", color="C0", lw=1.6, label="temperature")
    axc.set_xlabel("quench step (record index)"); axc.set_ylabel("normalised")
    axc.set_title("the fluid condenses: energy & baryon content collapse, a relic freezes", fontsize=9.5)
    axc.legend(fontsize=8, loc="upper right"); axc.grid(alpha=0.2)

    axt = fig.add_subplot(gs[1, 2:])
    axt.plot(t, res["B"], "o-", color="C2", ms=3)
    axt.axhline(0, color="0.6", lw=0.8, ls=":")
    axt.set_xlabel("quench step"); axt.set_ylabel("net baryon number $B$")
    axt.set_title("net baryon number (topologically conserved through the quench)", fontsize=9.5)
    axt.grid(alpha=0.2)

    fig.suptitle("3D genesis: baryons (torsitons) crystallising from the cooling chiral fluid "
                 "(SU(2) Skyrme quench)", fontsize=12, y=1.0)
    out = os.path.join(PDF_DIR, "genesis_quench3d.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  final: content={res['content'][-1]:.1f}  netB={res['B'][-1]:.2f}  "
          f"E {res['energy'][0]:.0f}->{res['energy'][-1]:.0f}")


if __name__ == "__main__":
    main()
