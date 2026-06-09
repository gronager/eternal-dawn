#!/usr/bin/env python3
r"""Horizon 3d: the unified genesis -- ALL observable matter condensing from ONE transition.

The capstone. As the chiral condensate v(T) (grey) switches on through a single transition, every
species acquires mass AT ONCE and the matter content of the universe fans out: the colourless leptons
(e, mu, tau) and the sub-eV neutrinos (3c), and the coloured sector confined (3a) into the baryon
octet and decuplet (3b). One condensate, two organising principles, the whole Standard-Model matter
spectrum -- 11 orders of magnitude -- out of one bounce. Renders figures/pdf/genesis_horizon.pdf.

The single owed input is the scale Lambda and the rung s_T: this runs on the calibrated values now and
RE-RUNS for the absolute-MeV prediction the moment the lattice gives them (genesis_horizon.genesis()).

HONEST: baryons ~30%; leptons from the soliton ladder + one s_T; neutrino values illustrative; the
structure (one transition, all species) is the content, the numbers are the lattice's word.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_horizon as gh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

SECTOR = {"octet": ("C0", "baryon octet"), "decuplet": ("C3", "baryon decuplet"),
          "lepton": ("C2", "charged leptons"), "neutrino": ("0.55", "neutrinos")}
SYM = {"N": "N", "Lambda": "Λ", "Sigma": "Σ", "Xi": "Ξ", "Delta": "Δ", "Sigma*": "Σ*",
       "Xi*": "Ξ*", "Omega": "Ω", "e": "e", "mu": "μ", "tau": "τ",
       "nu1": "ν₁", "nu2": "ν₂", "nu3": "ν₃"}


def main():
    g = gh.genesis()
    temps, v, masses, matter = g["temps"], g["v"], g["masses"], g["matter"]
    floor = 1e-8

    fig, ax = plt.subplots(figsize=(11.6, 6.8))
    for sp, mt in masses.items():
        sec = matter[sp][1]
        col = SECTOR[sec][0]
        ax.plot(temps, np.clip(mt, floor, None), color=col, lw=2.0,
                alpha=0.9 if sec != "neutrino" else 0.6)
        ax.annotate(SYM[sp], (temps[-1], max(matter[sp][0], floor)), xytext=(4, 0),
                    textcoords="offset points", color=col, fontsize=9, va="center", fontweight="bold")
    # the condensate
    axv = ax.twinx()
    axv.plot(temps, v, color="0.4", lw=1.6, ls=":")
    axv.set_ylabel("condensate $v(T)$", color="0.4"); axv.set_ylim(0, 1.3)
    axv.tick_params(axis="y", labelcolor="0.4")

    ax.axvline(g["Tc"], color="0.7", lw=1.0)
    ax.text(g["Tc"], floor * 2, " one transition\n (all matter switches on)", fontsize=8.5,
            color="0.4", ha="left", va="bottom")
    ax.set_yscale("log"); ax.set_xlim(1.0, 0.0); ax.set_ylim(floor, 1e5)
    ax.set_xlabel("temperature  $T$   (cooling $\\rightarrow$)")
    ax.set_ylabel("mass  (MeV)")
    ax.set_title("Horizon 3d — the unified genesis: ALL observable matter from one chiral transition\n"
                 "leptons (colourless rungs) + neutrinos (weakest grip) + baryons (coloured Skyrmions, "
                 "confined)", fontsize=10.8)
    handles = [plt.Line2D([], [], color=c, lw=2.5, label=lab) for c, lab in SECTOR.values()]
    handles.append(plt.Line2D([], [], color="0.4", lw=1.6, ls=":", label="condensate $v(T)$"))
    ax.legend(handles=handles, fontsize=8.5, loc="center left")
    ax.grid(alpha=0.16, which="both")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "genesis_horizon.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    span = max(m for m, _ in matter.values()) / min(m for m, _ in matter.values())
    print(f"  {len(matter)} species, mass span {span:.1e} ({np.log10(span):.0f} orders)")
    print("  RE-RUN hook ready: genesis_horizon.genesis(scale_MeV=Lambda, s_T=s_T) once the lattice lands")


if __name__ == "__main__":
    main()
