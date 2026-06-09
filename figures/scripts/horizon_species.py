#!/usr/bin/env python3
r"""The Horizon, assembled: the whole observable fermion spectrum from ONE chiral soliton.

Every observable matter particle is the same chiral soliton in a different gauge sector:
  * LEPTONS (colourless, free): the soliton's RADIAL rungs -> the generations e, mu, tau (3c);
    the neutrinos are the same rungs with the weakest condensate grip -> sub-eV (3c).
  * BARYONS (coloured, confined): the topological Skyrmion + FLAVOUR rotation -> the octet
    (N, Lambda, Sigma, Xi) and decuplet (Delta, Sigma*, Xi*, Omega) (3b).
  * The quarks themselves are CONFINED inside the baryons (3a) -- never free.
One condensate, one configurational mass, two organising principles set by colour. Plotted on a log
mass axis spanning the sub-eV neutrino to the GeV Omega -- the Standard Model's matter, condensing out
of the bounce. Renders figures/pdf/horizon_species.pdf.

HONEST: baryons ~30% (Skyrme + rigid rotor, 3b); charged leptons from the soliton ladder + one s_T
(3c); neutrino masses illustrative (the mechanism, not the values); absolute scale calibrated, becomes
first-principles with Lambda from the lattice.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_su3 as su3
from cartasis_sims import genesis_leptons as lep

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    bar = su3.baryon_spectrum()["mass"]                     # MeV
    L = lep.lepton_spectrum()
    charged = dict(zip(["e", "μ", "τ"], L["charged"]["mass"]))
    nu = dict(zip(["ν₁", "ν₂", "ν₃"], L["neutrino"]["mass"]))

    # columns: (x position, colour, label, {name: mass MeV}, organising principle)
    cols = [
        (0, "0.55", "neutrinos\n(weakest grip)", nu, "radial rungs"),
        (1, "C2", "charged leptons\n$e,\\mu,\\tau$", charged, "radial rungs (generations)"),
        (2, "C0", "baryon octet\n$J=1/2$", {"N": bar["N"], "Λ": bar["Lambda"],
                                            "Σ": bar["Sigma"], "Ξ": bar["Xi"]}, "Skyrmion + flavour"),
        (3, "C3", "baryon decuplet\n$J=3/2$", {"Δ": bar["Delta"], "Σ*": bar["Sigma*"],
                                               "Ξ*": bar["Xi*"], "Ω": bar["Omega"]}, "Skyrmion + flavour"),
    ]

    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    for xpos, col, lab, group, _ in cols:
        for nm, mass in group.items():
            ax.plot([xpos - 0.28, xpos + 0.28], [mass, mass], "-", color=col, lw=2.4)
            ax.annotate(nm, (xpos + 0.30, mass), fontsize=10, color=col, va="center", fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks([c[0] for c in cols]); ax.set_xticklabels([c[1 + 1] for c in cols], fontsize=9.5)
    for xpos, col, lab, group, _ in cols:
        ax.text(xpos, 4e3, lab, ha="center", va="bottom", fontsize=9.5, color=col)
    ax.set_xlim(-0.6, 3.7); ax.set_ylim(1e-8, 1e5)
    ax.set_ylabel("mass  (MeV)")
    ax.set_xticks([])
    # the two organising principles
    ax.axvspan(-0.55, 1.55, color="C2", alpha=0.045)
    ax.axvspan(1.55, 3.6, color="C0", alpha=0.045)
    ax.text(0.5, 2e-7, "COLOURLESS · free\nradial rungs $\\to$ generations", ha="center", fontsize=9,
            color="C2", style="italic")
    ax.text(2.5, 2e-7, "COLOURED · confined (3a)\nSkyrmion + flavour $\\to$ octet", ha="center",
            fontsize=9, color="C0", style="italic")
    ax.set_title("The observable fermion spectrum from ONE chiral soliton\n"
                 "leptons = colourless rungs (3c)   ·   baryons = coloured Skyrmions (3b)   ·   "
                 "quarks confined within (3a)", fontsize=11.5)
    ax.grid(alpha=0.18, axis="y", which="both")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "horizon_species.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  span: neutrino {min(nu.values()):.2e} MeV  ->  Omega {bar['Omega']:.0f} MeV "
          f"({np.log10(bar['Omega'] / min(nu.values())):.0f} orders of magnitude)")


if __name__ == "__main__":
    main()
