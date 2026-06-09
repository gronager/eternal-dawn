#!/usr/bin/env python3
r"""Genesis at the particle scale: torsitons condensing from a cooling chiral fluid (2D prototype).

Top row: snapshots of the topological-charge density q(x) as the hot, disordered post-bounce "tally"
is cooled through the chiral transition -- the torsitons (spinning winding-knots) crystallise out as
sharp +/- peaks. Bottom-left: the condensation -- energy and the unsigned winding content sum|q| (the
torsiton-count proxy) falling as the fluid orders and +/- pairs annihilate, leaving a frozen relic
set. Bottom-right: the Kibble--Zurek signature -- a FASTER quench freezes in MORE torsitons (the
correlation length can't keep up), the same universal mechanism behind cosmic strings and superfluid
vortices. Renders figures/pdf/genesis_quench.pdf.

HONEST: classical stochastic effective field theory (the genesis MECHANISM + defect density + KZ
scaling), not the exact quantum masses; 2D baby-Skyrmions, dimensionless units. The picture, runnable.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_quench as gq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

PAR = dict(kappa=5.0, mu2=0.01, dt=0.03, gamma=0.5)         # stable resolved baby-Skyrmions (size ~5)


def main():
    res = gq.run_quench(L=160, steps=6000, T0=0.9, T1=0.0, seed=7, record=40, **PAR)
    frames = res["frames"]
    # pick four snapshots: early (disordered), two intermediate, final (frozen torsitons)
    idx = [int(0.08 * len(frames)), int(0.30 * len(frames)), int(0.60 * len(frames)), len(frames) - 1]

    fig = plt.figure(figsize=(12.6, 8.4))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.12, 1.12, 1.0], hspace=0.34, wspace=0.16)

    qmax = np.percentile(np.abs(frames[idx[1]][2]), 99.7) or 1.0
    labels = ["hot tally\n(disordered)", "cooling\n(domains)", "condensing\n(knots forming)",
              "frozen relic\n(torsitons)"]
    for k, i in enumerate(idx):
        s, n, q = frames[i]
        T = res["T"][i]
        # row 0: the CONDENSATE itself -- the medium ordering from disordered (~0) to vacuum (~1)
        axo = fig.add_subplot(gs[0, k])
        axo.imshow(gq.coarse_order(n, sigma=3.0).T, origin="lower", cmap="viridis", vmin=0, vmax=1,
                   interpolation="nearest")
        axo.set_xticks([]); axo.set_yticks([])
        axo.set_title(f"{labels[k]}\nT={T:.2f}", fontsize=9.5)
        if k == 0:
            axo.set_ylabel("condensate $|\\langle n\\rangle|$", fontsize=9.5)
        # row 1: the torsitons -- the topological charge density (the knots in that medium)
        axq = fig.add_subplot(gs[1, k])
        axq.imshow(q.T, origin="lower", cmap="RdBu_r", vmin=-qmax, vmax=qmax, interpolation="nearest")
        axq.set_xticks([]); axq.set_yticks([])
        if k == 0:
            axq.set_ylabel("torsitons $q(x)$", fontsize=9.5)

    # --- condensation: energy + winding content vs time ---
    axc = fig.add_subplot(gs[2, :2])
    t = np.arange(len(res["content"]))
    axc.plot(t, res["content"] / res["content"].max(), "o-", color="C3", ms=3,
             label="winding content $\\sum|q|$ (torsiton count)")
    axc.plot(t, res["energy"] / res["energy"].max(), "s-", color="0.5", ms=3, label="energy")
    axc.plot(t, res["T"] / res["T"].max(), ":", color="C0", lw=1.6, label="temperature")
    axc.set_xlabel("quench step (record index)"); axc.set_ylabel("normalised")
    axc.set_title("the fluid condenses: energy & winding fall, a relic set freezes", fontsize=10)
    axc.legend(fontsize=8, loc="upper right"); axc.grid(alpha=0.2)

    # --- Kibble-Zurek: faster quench -> more torsitons ---
    axk = fig.add_subplot(gs[2, 2:])
    rates, content = gq.quench_rate_scan([0.5, 1.0, 2.0, 4.0, 8.0], L=96, base_steps=8000,
                                         T0=0.9, seed=11, **PAR)
    axk.loglog(rates, content, "o-", color="C2", lw=1.8)
    axk.set_xlabel("cooling rate  (faster $\\rightarrow$)"); axk.set_ylabel("torsitons frozen in")
    axk.set_title("Kibble--Zurek: faster quench, more torsitons", fontsize=10)
    axk.grid(alpha=0.2, which="both")

    fig.suptitle("Particle-scale genesis: torsitons crystallising from the cooling chiral fluid "
                 "(2D baby-Skyrme quench)", fontsize=12, y=0.99)
    out = os.path.join(PDF_DIR, "genesis_quench.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  final: content={res['content'][-1]:.1f}  netQ={res['Q'][-1]:.1f}  "
          f"E {res['energy'][0]:.0f}->{res['energy'][-1]:.0f}")
    print(f"  KZ rates {rates} -> content {np.round(content,1)}")


if __name__ == "__main__":
    main()
