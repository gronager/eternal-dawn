#!/usr/bin/env python3
r"""The Horizon: the twelve fermions condensing out of ONE chiral transition.

As the universe cools and the condensate v(T) switches on (the same S-curve the genesis quench
produces), every fermion acquires mass AT THE SAME transition -- the framework's unification, unlike
the Standard Model's staged electroweak-then-QCD -- and fans out into the observed hierarchy. Each
mass factorises as m = Lambda * c_T * |O_n(s_T)|^2: the tower coupling c_T (colour/charge) sets where
a tower sits, the radial rung O_n(s_T) sets the three generations. Coloured by tower; the condensate
curve is the grey S. Renders figures/pdf/horizon_spectrum.pdf.

HONEST: RELATIVE masses (absolute MeV needs Lambda = g,b0); the inter-tower factors and s_T are the
lattice residual -- the framework's content is that the spread FACTORISES as c_T x rung at one scale.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import horizon as hz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

COL = {"up-quark": "C3", "down-quark": "C0", "charged-lepton": "C2", "neutrino": "0.6"}
NAME = {"up-quark": ["u", "c", "t"], "down-quark": ["d", "s", "b"],
        "charged-lepton": ["e", "μ", "τ"], "neutrino": ["ν₁", "ν₂", "ν₃"]}


def main():
    temps, v, masses, rel = hz.emergence(n_T=300, Tc=0.5, width=0.05)
    floor = 1e-2                                              # log-plot floor (relative units)

    fig, ax = plt.subplots(figsize=(11.0, 6.4))
    for t in hz.TOWERS:
        m = np.clip(masses[t], floor, None)
        for g in range(3):
            ax.plot(temps, m[:, g], color=COL[t], lw=2.0,
                    alpha=0.9 if t != "neutrino" else 0.55,
                    ls="-" if t != "neutrino" else "--")
            ax.annotate(NAME[t][g], (temps[-1], m[-1, g]), xytext=(4, 0),
                        textcoords="offset points", color=COL[t], fontsize=8.5, va="center")
    # the condensate S-curve (right axis, the medium switching on)
    axv = ax.twinx()
    axv.plot(temps, v, color="0.4", lw=1.6, ls=":")
    axv.fill_between(temps, 0, v, color="0.5", alpha=0.06)
    axv.set_ylabel("condensate  $v(T)$", color="0.4"); axv.set_ylim(0, 1.25)
    axv.tick_params(axis="y", labelcolor="0.4")

    ax.axvline(0.5, color="0.7", lw=1.0, ls="-")
    ax.text(0.5, floor * 1.4, "one transition\n(all masses switch on)", fontsize=8.5, color="0.4",
            ha="center", va="bottom")
    ax.set_yscale("log")
    ax.set_xlim(temps[0], temps[-1])                        # hot (T=1) left -> cold (T=0) right
    ax.set_xlabel("temperature  $T$   (cooling $\\rightarrow$)")
    ax.set_ylabel("fermion mass   (relative, electron $=1$)")
    ax.set_title("The Horizon: the twelve fermions condensing out of one chiral transition\n"
                 "$m = \\Lambda\\, c_T\\, |O_n(s_T)|^2$ — tower coupling $c_T$ (colour) sets the bands, "
                 "rung $s_T$ sets the generations", fontsize=10.5)
    # tower legend
    handles = [plt.Line2D([], [], color=COL[t], lw=2.5,
                          label={"up-quark": "up-type (u,c,t)", "down-quark": "down-type (d,s,b)",
                                 "charged-lepton": "leptons (e,μ,τ)",
                                 "neutrino": "neutrinos (illustrative)"}[t]) for t in hz.TOWERS]
    ax.legend(handles=handles, fontsize=8.5, loc="lower left")
    ax.grid(alpha=0.18, which="both")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "horizon_spectrum.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    for t in hz.TOWERS:
        print(f"  {t:14s}: c_T={rel[t]['c_T']:.2e}  s_T={rel[t]['s_T']:.3f}  "
              f"final masses(rel)={np.round(rel[t]['m'], 2)}")


if __name__ == "__main__":
    main()
