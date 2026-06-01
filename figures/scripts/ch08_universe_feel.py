#!/usr/bin/env python3
r"""What other universe types are like: baryon-rich BHUs, and observer-free OGUs.

Left: the accretion-family axis -- our clean, photon-bathed cosmos (eta ~ 1 ppb) vs a
baryon-rich/degenerate one (eta -> 1): a short-lived furnace with no radiation era.
Right: structure growth with vs without dark matter -- an OGU has no parent, hence no
projected dark matter, hence baryon-only growth that falls ~100x short of forming
galaxies. So OGUs make black holes but essentially no observers: the census shows none
by prediction, not omission. Renders figures/pdf/universe_feel.pdf, writes
sims/output/universe_feel.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import universe_feel as uf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    short = (uf.delta_today(uf.growth_factor_baryon_only(), uf.SEED_BARYON_ONLY)
             / uf._COLLAPSE_BAR)

    lines = [
        "Universe types: the baryon-rich feel, and why OGUs have no observers",
        "=" * 68,
        "A. BARYON-RICH BHU (eta -> 1, fed by concentrated baryons):",
        "   - photons per baryon ~ 1 (vs ~1.6e9 for us): NO radiation era, no CMB;",
        "   - gas cools and fragments at once -> one fast collapse into compact",
        "     objects/black holes, not a slow cosmic web;",
        "   - a degenerate, short-lived furnace -- almost certainly sterile.",
        "   Our clean, photon-bathed, slow cosmos is the opposite extreme.",
        "",
        "B. OGU (forms from the void, NO parent):",
        "   - no parent -> no projected mass -> NO dark matter, no inherited lumpiness;",
        "   - baryon-only structure grows only from recombination, from a Silk-damped",
        "     seed: it clears only ~{:.0f}% of our universe's collapse bar;".format(short*100),
        f"   - forms galaxies? {uf.forms_galaxies(False)}.  hosts observers? "
        f"{uf.ogu_hosts_observers()}.",
        "   An OGU is good at exactly ONE thing -- making the rare black hole (hence",
        "   descendants) -- and almost nothing else: smooth, homogeneous, observer-free.",
        "",
        "READING: this is a SECOND, independent reason we are not an OGU. We already",
        "knew it from observing dark matter + dark energy (which an OGU lacks); now we",
        "see an OGU could not host us anyway -- with no dark matter it stays nearly",
        "homogeneous and builds no galaxies. So the census correctly shows NO OGU",
        "observers: that is a prediction, not an omission. Observers live in BHUs,",
        "in the clean family -- exactly our footprint.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "universe_feel.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: photons per baryon across the family axis.
    etas = np.logspace(-10, 0, 200)
    axL.plot(etas, 1.0 / etas, "C0", lw=2.2)
    axL.axvline(uf.ETA_US, color="C2", lw=1.6)
    axL.text(uf.ETA_US * 1.3, 1e3, "us (clean):\n$\\eta\\sim$1 ppb,\nphoton-bathed,\nslow",
             fontsize=8.5, color="C2")
    axL.axvspan(1e-3, 1.0, color="C3", alpha=0.12)
    axL.text(3e-3, 3e6, "baryon-rich / degenerate:\nno radiation era,\nfast collapse,\nsterile",
             fontsize=8.5, color="C3")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"baryon-to-photon ratio  $\eta$  (accretion family)")
    axL.set_ylabel(r"photons per baryon  $1/\eta$")
    axL.set_title("What a baryon-rich BHU feels like\n"
                  "lose the photon bath $\\Rightarrow$ degenerate furnace", fontsize=11)
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: structure growth, OGU vs BHU.
    cats = ["OGU\n(no dark matter)", "BHU\n(with dark matter,\nlike us)"]
    vals = [uf.delta_today(uf.growth_factor_baryon_only(), uf.SEED_BARYON_ONLY),
            uf._COLLAPSE_BAR]
    vals = [v / uf._COLLAPSE_BAR for v in vals]      # normalize to our bar = 1
    colors = ["0.6", "C0"]
    axR.bar(cats, vals, color=colors, alpha=0.85, width=0.6)
    axR.axhline(1.0, color="C3", ls="--", lw=1.4)
    axR.text(0.5, 1.1, "collapse bar (galaxies form)", fontsize=8.5, color="C3",
             ha="center")
    for c, v in zip(cats, vals):
        axR.text(c, v + 0.04, f"{v:.2f}", ha="center", fontsize=9)
    axR.text(0, 0.4, "structure-poor\nNO galaxies\nNO observers", ha="center",
             fontsize=8.5, color="0.3")
    axR.set_ylabel(r"structure growth  (our universe $=1$)")
    axR.set_title("Why OGUs have no observers\n"
                  "no parent $\\Rightarrow$ no dark matter $\\Rightarrow$ no galaxies",
                  fontsize=11)
    axR.set_ylim(0, 1.3)
    axR.grid(True, axis="y", alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "universe_feel.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "universe_feel.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'universe_feel.pdf')}")


if __name__ == "__main__":
    main()
