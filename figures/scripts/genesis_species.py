#!/usr/bin/env python3
r"""Species-resolved genesis (Horizon step 2): distinct species condensing from ONE transition, masses
from the couplings. Several chiral SECTORS (different charge/colour couplings) are batched and cooled
together. Left: the species rest mass from a relaxed single knot -- M ~ sqrt(kappa), the hierarchy
coming OUT of the soliton dynamics given the couplings. Middle: all sectors condense at ONE transition
(the unification), their relic baryon content fanning by coupling. Right: the final baryon density of
the lightest vs heaviest sector -- distinct knots, distinct sizes. Renders genesis_species.pdf.

HONEST: a schematic multi-sector model (each tower an O(4) Skyrme sector, coupling from c_T) over a
MODEST spread -- the full x3477 hierarchy needs the s_T structure (the lattice residual), too extreme
for one classical box. It shows the mechanism: distinct species, masses from couplings, one transition.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_species as gs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

SPECIES = [("sector 1", 12.0, 0.015), ("sector 2", 20.0, 0.011),
           ("sector 3", 32.0, 0.008), ("sector 4", 48.0, 0.006)]
COLS = ["C2", "C0", "C1", "C3"]


def main():
    # (1) species rest masses from relaxed single knots: M ~ sqrt(kappa)
    relax = gs.relax_solitons(SPECIES, L=40, steps=600, dt=0.004)
    # (2) the condensation: all sectors cool through one transition, batched
    run = gs.run_species(SPECIES, L=44, steps=1600, dt=0.005, seed=4, record=16)

    fig, (axM, axC, axB) = plt.subplots(1, 3, figsize=(13.4, 4.4),
                                        gridspec_kw={"width_ratios": [1.0, 1.1, 0.9]})

    # --- left: mass from coupling (M ~ sqrt(kappa)) ---
    kap = relax["kappa"]
    ks = np.linspace(kap.min(), kap.max(), 50)
    axM.plot(ks, relax["mass"][0] * np.sqrt(ks / kap[0]), "k--", lw=1.2, alpha=0.6,
             label=r"$M\propto\sqrt{\kappa}$ (Skyrme)")
    for k, m, c, nm in zip(kap, relax["mass"], COLS, relax["names"]):
        axM.plot(k, m, "o", color=c, ms=10)
        axM.annotate(nm, (k, m), xytext=(6, -2), textcoords="offset points", fontsize=8.5, color=c)
    axM.set_xlabel(r"coupling $\kappa$  (from $c_T$)"); axM.set_ylabel("species rest mass (knot energy)")
    axM.set_title("masses OUT of the dynamics:\nheavier coupling, heavier knot", fontsize=10)
    axM.legend(fontsize=8); axM.grid(alpha=0.2)

    # --- middle: condensation, all at one transition ---
    t = np.arange(run["content"].shape[0])
    for s in range(len(SPECIES)):
        axC.plot(t, run["content"][:, s] / run["content"][:, s].max(), "o-", color=COLS[s], ms=3,
                 label=SPECIES[s][0])
    axC.plot(t, run["T"] / run["T"].max(), ":", color="0.5", lw=1.6, label="temperature")
    axC.set_xlabel("quench step (record)"); axC.set_ylabel("baryon content (norm.)")
    axC.set_title("all species condense at ONE transition\n(the unification)", fontsize=10)
    axC.legend(fontsize=7.5, loc="upper right"); axC.grid(alpha=0.2)

    # --- right: distinct knots, lightest vs heaviest sector (final baryon density projection) ---
    b = run["final_b"]                                       # (B, L, L, L)
    light = np.max(np.abs(b[0]), axis=2)
    heavy = np.max(np.abs(b[-1]), axis=2)
    vmax = np.percentile(np.concatenate([light.ravel(), heavy.ravel()]), 99.7) or 1.0
    axB.imshow(np.concatenate([light, np.full((light.shape[0], 2), np.nan), heavy], axis=1).T,
               origin="lower", cmap="inferno", vmin=0, vmax=vmax, interpolation="nearest")
    axB.set_xticks([]); axB.set_yticks([])
    axB.set_title(f"distinct knots:\n{SPECIES[0][0]} (left) vs {SPECIES[-1][0]} (right)", fontsize=10)

    fig.suptitle("Species-resolved genesis: distinct species condensing from one transition, masses "
                 "from the couplings  (batched Skyrme, GPU-ready)", fontsize=11.5, y=1.02)
    fig.tight_layout()
    out = os.path.join(PDF_DIR, "genesis_species.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}   (backend GPU={run['gpu']})")
    print(f"  species masses (rel): {np.round(relax['mass_rel'], 2)}  for kappa {relax['kappa']}")


if __name__ == "__main__":
    main()
