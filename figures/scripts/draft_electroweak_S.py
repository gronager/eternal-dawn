#!/usr/bin/env python3
r"""The electroweak S parameter from the soliton sector -- leading order, honestly.

A: the constituent-fermion vector and axial spectral functions and their chiral-odd
   difference rho_V - rho_A > 0 -- why S > 0 (the vector channel outweighs the axial).
B: the verdict -- the leading-order S (constituent loop, S = N_c/6pi per doublet) sits
   ABOVE the LEP bound (the graveyard) for QCD-like and for Pati-Salam (N_c=4) content;
   escaping needs the resonances to 'walk' (M_V/f_pi from ~8 to >~14), the subleading
   correction this calculation does not contain.

Renders figures/pdf/electroweak_S.pdf, writes sims/output/electroweak_S.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import electroweak_S as ew

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    M = 1.0
    contents = [("QCD-like ($N_c$=3, 1 doublet)", 3, 1.0, "C3"),
                ("2 doublets ($N_c$=3)", 3, 2.0, "C1"),
                ("Pati-Salam ($N_c$=4)", 4, 1.0, "C4")]

    lines = ["The electroweak S parameter from the soliton sector (leading order)",
             "=" * 66,
             "  constituent-loop S = N_c n_doublets / (6 pi)  (Peskin-Takeuchi):"]
    for name, nc, nd, _ in contents:
        v = ew.verdict(nc, nd)
        lines.append(f"    {name:32s} S_leading = {v['S_leading']:.3f}  "
                     f"(> LEP bound {ew.S_LEP_BOUND}: graveyard)")
    v = ew.verdict(3, 1)
    lines += [
        "",
        f"  to escape (S < 0.1): need M_V/f_pi > {v['MV_over_fpi_needed']:.1f}  "
        f"(vs QCD's {v['qcd_MV_over_fpi']:.1f}) -- a ~1.7x 'walk'.",
        "",
        "VERDICT (honest): the leading-order S from our model parameters is in the",
        "graveyard (~0.16-0.32), and Pati-Salam's extra colour makes it WORSE, not",
        "better. The escape needs subleading resonance ('walking') corrections -- the V",
        "and A composites rearranging so their spectral functions cancel -- which need",
        "the RPA correlators this calculation does not contain. So S stays the owed",
        "make-or-break; degeneracy/Pauli says only that the DIRECTION is favourable.",
        "Part III can still die here, exactly as its own status admits.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "electroweak_S.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0))

    # ---- A: spectral functions ----
    s = np.linspace(4 * M**2, 60 * M**2, 800)
    axA.plot(s, ew.spectral_V(s, M), "C0", lw=2.0, label=r"$\rho_V(s)$ (vector)")
    axA.plot(s, ew.spectral_A(s, M), "C1", lw=2.0, label=r"$\rho_A(s)$ (axial)")
    axA.fill_between(s, ew.spectral_difference(s, M), color="C3", alpha=0.25)
    axA.plot(s, ew.spectral_difference(s, M), "C3", lw=2.0,
             label=r"$\rho_V-\rho_A>0$  (drives $S>0$)")
    axA.set_xlabel(r"$s/M^2$")
    axA.set_ylabel(r"spectral function")
    axA.set_title("Constituent loop: vector outweighs axial\n"
                  r"$\Rightarrow S>0$ (the graveyard sign)", fontsize=11)
    axA.legend(fontsize=9, loc="upper right")
    axA.grid(True, alpha=0.2)
    axA.set_xlim(0, 60)

    # ---- B: the verdict landscape ----
    MV_over_fpi = np.linspace(5, 20, 200)
    S = np.array([ew.s_walking(1.0 / x) for x in MV_over_fpi])
    axB.plot(MV_over_fpi, S, "C2", lw=2.0, label="walking curve (Weinberg sum rules)")
    axB.axhline(ew.S_LEP_BOUND, color="0.4", ls="--", lw=1.2)
    axB.text(5.2, 0.11, "LEP bound", fontsize=8.5, color="0.35")
    axB.axhspan(ew.S_LEP_BOUND, 0.45, color="0.85", alpha=0.6)
    for name, nc, nd, c in contents:
        S0 = ew.s_leading(nc, nd)
        axB.plot(ew.verdict(nc, nd)["qcd_MV_over_fpi"], S0, "o", color=c, ms=9, label=name)
    need = ew.verdict(3, 1)["MV_over_fpi_needed"]
    axB.plot(need, ew.S_LEP_BOUND, "kv", ms=9)
    axB.annotate("escape needs\n'walking'", (need, ew.S_LEP_BOUND),
                 textcoords="offset points", xytext=(-10, 18), fontsize=8.5, ha="center")
    axB.set_xlim(5, 20); axB.set_ylim(0, 0.45)
    axB.set_xlabel(r"$M_V/f_\pi$  (resonance heaviness / walking)")
    axB.set_ylabel("Peskin-Takeuchi $S$")
    axB.set_title("Leading $S$ is in the graveyard;\nescape needs subleading walking",
                  fontsize=11)
    axB.legend(fontsize=8, loc="upper right")
    axB.grid(True, alpha=0.2)

    fig.suptitle("The electroweak S: leading order from the model is the graveyard; "
                 "the escape (walking) is the owed make-or-break", fontsize=12, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "electroweak_S.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "electroweak_S.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'electroweak_S.pdf')}")


if __name__ == "__main__":
    main()
