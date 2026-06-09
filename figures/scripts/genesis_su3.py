#!/usr/bin/env python3
r"""Horizon 3b: the Skyrme baryon spectrum -- the nucleon, the Delta, and the hyperons from ONE
soliton. Left: the hedgehog chiral angle F(r) (the soliton profile, F(0)=pi -> 0). Right: the baryon
octet (J=1/2) and decuplet (J=3/2) vs strangeness -- predicted (filled, lines) against observed (open
stars). N, Delta, Omega calibrate (f_pi, e, the strange scale); the rest are PREDICTIONS, landing
within ~1% (decuplet equal spacing) to ~9% (octet). The first piece of the Horizon that is a
prediction, not a picture. Renders figures/pdf/genesis_su3.pdf.

HONEST: classical soliton + rigid-rotor collective quantization (~30% Skyrme accuracy). The octet
Lambda-Sigma split (their different isospin) needs the full SU(3)/Yabu-Ando treatment -- here they are
degenerate. With f_pi,e from the lattice Lambda this becomes absolute-scale from first principles.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis_su3 as su3

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    x, F, M_red = su3.hedgehog_profile()
    res = su3.baryon_spectrum()
    m, obs, s = res["mass"], res["obs"], res["scales"]

    fig, (axF, axS) = plt.subplots(1, 2, figsize=(12.0, 4.8), gridspec_kw={"width_ratios": [1.0, 1.5]})

    # --- left: the hedgehog profile ---
    axF.plot(x, F, "C0-", lw=2.2)
    axF.axhline(np.pi, color="0.7", ls=":", lw=1); axF.axhline(0, color="0.7", ls=":", lw=1)
    axF.set_xlim(0, 8); axF.set_ylim(-0.1, np.pi + 0.2)
    axF.set_yticks([0, np.pi / 2, np.pi]); axF.set_yticklabels(["0", "π/2", "π"])
    axF.set_xlabel("radius  $r$  (soliton units)"); axF.set_ylabel("chiral angle  $F(r)$")
    axF.set_title(f"the baryon soliton (hedgehog)\n$M_{{cl}}$={s['M_cl']:.0f} MeV, "
                  f"$f_\\pi$={s['f_pi']:.0f}, $e$={s['e']:.2f}", fontsize=10)
    axF.grid(alpha=0.2)

    # --- right: the octet + decuplet spectrum vs strangeness ---
    octet = ["N", "Lambda", "Sigma", "Xi"]
    decup = ["Delta", "Sigma*", "Xi*", "Omega"]
    syms = {"N": "N", "Lambda": "Λ", "Sigma": "Σ", "Xi": "Ξ",
            "Delta": "Δ", "Sigma*": "Σ*", "Xi*": "Ξ*", "Omega": "Ω"}
    for group, col, lab, J in [(octet, "C2", "octet  $J=1/2$", "1/2"), (decup, "C3", "decuplet  $J=3/2$", "3/2")]:
        sp = [su3.STRANGENESS[b] for b in group]
        pm = [m[b] for b in group]
        om = [obs[b] for b in group]
        axS.plot(sp, pm, "o-", color=col, lw=1.8, ms=9, label=f"{lab} (predicted)")
        axS.plot(sp, om, "*", color=col, ms=15, mfc="white", mec=col, mew=1.6,
                 label=f"{lab} (observed)")
        for b, ss, p in zip(group, sp, pm):
            cal = b in ("N", "Delta", "Omega")
            axS.annotate(syms[b] + ("*" if cal else ""), (ss, p), xytext=(7, -3),
                         textcoords="offset points", color=col, fontsize=10, fontweight="bold")
    axS.set_xlabel("strangeness  $|S|$"); axS.set_ylabel("baryon mass  (MeV)")
    axS.set_xticks([0, 1, 2, 3])
    axS.set_title("baryon octet & decuplet: predicted vs observed\n"
                  "(* = calibration: N, Δ, Ω; the rest predicted within ~1–9%)", fontsize=10)
    axS.legend(fontsize=8, loc="upper left"); axS.grid(alpha=0.2)

    fig.suptitle("Horizon 3b: the baryon spectrum from one Skyrme soliton — the Horizon becomes a "
                 "prediction", fontsize=11.5, y=1.0)
    fig.tight_layout()
    out = os.path.join(PDF_DIR, "genesis_su3.pdf")
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=130,
                                                       bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    for b in ["N", "Lambda", "Sigma", "Xi", "Delta", "Sigma*", "Xi*", "Omega"]:
        print(f"  {b:7s}: {m[b]:7.1f} vs {obs[b]:7.1f}  ({100*(m[b]-obs[b])/obs[b]:+.0f}%)")


if __name__ == "__main__":
    main()
