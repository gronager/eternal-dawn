#!/usr/bin/env python3
r"""Generations from the self-consistent torsiton well: arithmetic levels -> geometric masses.

Left: the lowest three radial rungs of the SELF-CONSISTENT bag (the generations) -- level energies
nearly evenly spaced (arithmetic), but configurational masses (the overlap of each rung with the
mass-giving substrate, Eq. configmass) that span orders of magnitude when the substrate is the SHARP
chiral-restored core (geometric) rather than the broad condensate (no hierarchy). The observed
charged leptons are overlaid. Right: the sharp-overlap spread grows with the well depth (coupling g);
the self-consistent wells (no fit) reach ~50, with the remaining climb to the observed lepton spread
(3477) the deeper chiral dynamics owed to the lattice. Renders figures/pdf/generations_overlap.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import generations as gen

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

LEPTONS = np.array([0.511, 105.66, 1776.86])     # e, mu, tau (MeV)


def main():
    h = gen.self_consistent_hierarchy(g=14.0, m_sigma=0.5)
    lv = h["levels"]
    broad = h["broad_masses"]
    sharp = h["sharp_masses"]

    # depth trend: sharp spread vs coupling g (converged 3-rung wells)
    gs, spreads = [], []
    for g in (10.0, 12.0, 14.0, 16.0):
        hh = gen.self_consistent_hierarchy(g=g, m_sigma=0.5)
        if hh["n_bound"] >= 3 and hh["converged"]:
            gs.append(g); spreads.append(hh["sharp_spread"])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.8, 4.0))

    idx = np.arange(3)
    axL.semilogy(idx, lv / lv.min(), "s--", color="0.6", label="level energy $E_n$ (arithmetic)")
    axL.semilogy(idx, broad / broad.min(), "o-", color="C0", label=r"broad overlap (spread $\sim1$)")
    axL.semilogy(idx, sharp / sharp.min(), "o-", color="C3",
                 label=rf"sharp-core overlap (spread ${h['sharp_spread']:.0f}$)")
    axL.semilogy(idx, LEPTONS / LEPTONS.min(), "*-", color="C2", ms=11,
                 label="observed $e,\\mu,\\tau$ (3477)")
    axL.set_xticks(idx); axL.set_xticklabels(["gen I", "gen II", "gen III"])
    axL.set_ylabel("mass / lightest")
    axL.set_title("arithmetic levels $\\to$ geometric masses", fontsize=10.5)
    axL.legend(fontsize=7.8, loc="upper left")
    axL.grid(True, which="both", alpha=0.2)

    axR.plot(gs, spreads, "o-", color="C3", lw=1.8, label="self-consistent wells (no fit)")
    axR.axhline(3477, color="C2", ls="--", lw=1.2, label="observed leptons (3477)")
    axR.axhspan(20, 60, color="C3", alpha=0.08)
    axR.set_yscale("log")
    axR.set_xlabel(r"coupling $g$  (well depth $\to$)")
    axR.set_ylabel("sharp-overlap spread (gen III / gen I)")
    axR.set_title("the hierarchy steepens with depth", fontsize=10.5)
    axR.text(10.2, 70, "physical (self-consistent)\nwells reach $\\sim$50", fontsize=8, color="C3")
    axR.legend(fontsize=8, loc="center right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "generations_overlap.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  levels {np.round(lv,3)} (spacings {np.round(np.diff(lv),3)})")
    print(f"  broad spread {h['broad_spread']:.2f}  sharp spread {h['sharp_spread']:.1f}")
    print(f"  sharp masses (light=1) {np.round(sharp,1)}  ratios {np.round(sharp[1:]/sharp[:-1],1)}")
    print(f"  observed leptons ratios [206.8, 16.8], spread 3477")


if __name__ == "__main__":
    main()
