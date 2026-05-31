#!/usr/bin/env python3
r"""Growth of a universe-hole, and the generation depth it allows.

dM/dt = (P_H/c^2)[(M/M_crit)^4 - 1]: evaporate below M_crit, runaway (~M^2, g=2)
above. Depletion of the local patch freezes the final mass. The growth time
tau_gen ~ tau0(T_bath) caps the depth a supraverse of given age can reach,
D_max = T_supra/tau_gen -- closing the population loop to P(BHU_n).
Renders figures/pdf/growth.pdf, writes sims/output/growth.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import growth as g
from cartasis_sims import population as pop

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

HUBBLE_TIME_S = 4.35e17     # ~ present age of our universe, an illustrative unit


def main() -> None:
    ages = {"T_supra = 1 Hubble time": HUBBLE_TIME_S,
            "T_supra = 1e6 Hubble times": 1e6 * HUBBLE_TIME_S}

    lines = [
        "Growth of a universe-hole, and the depth it allows",
        "=" * 58,
        "  dM/dt = (P_H/c^2)[(M/M_crit)^4 - 1];  M_crit = hbar c^3/(8 pi G kB T_bath)",
        "  evaporate (M<M_crit) | runaway ~M^2 (g=2) above | depletion freezes it",
        "",
        f"  runaway is fast: t'(x0=2) = {g.runaway_time_universal(2):.3f} tau0,"
        f"  t'(x0=10) = {g.runaway_time_universal(10):.3f} tau0",
        f"  depletion freeze-out: M_final ~ x0 + bath (eats its causal patch)",
        "",
        "  void bath sets everything (tau_gen ~ 1/T_bath^3):",
        "    T_bath[K]    M_crit[kg]    tau_gen[s]",
    ]
    for T in (1e9, 1e12, 1e15):
        lines.append(f"    {T:.0e}    {g.m_crit(T):.2e}    {g.generation_time(T):.2e}")

    lines += ["", "  CLOSING THE LOOP -- depth cap D_max = T_supra/tau_gen, then P(BHU_n):"]
    for label, age in ages.items():
        lines.append(f"  {label}:")
        lines.append("    T_bath[K]   D_max   regime / P(BHU1 or 2)")
        for T in (1e11, 1e12, 1e13, 1e15):
            D = g.max_depth(age, T)
            if D == 0:
                note = "sterile foam (no BHUs) -- EXCLUDED by our existence"
            elif D <= 2:
                note = f"forced shallow: P(BHU1 or 2) = 1.00 (D_max={D})"
            else:
                # depth set by branching; quote a supercritical example m=1.8
                note = (f"depth set by m; e.g. m=1.8 -> "
                        f"P(BHU1 or 2)={pop.shallow_probability(1.8):.1e}")
            lines.append(f"    {T:.0e}     {D:<6d}  {note}")
    lines += [
        "",
        "READING: the void bath temperature is the master dial (tau_gen ~ 1/T^3).",
        "Too cold (T<~1e12 K here) and nothing grows in a Hubble-aged supraverse --",
        "a sterile foam of OGUs, excluded because we ARE a BHU (n>=1). A Goldilocks",
        "T_bath ~ 1e12 K gives D_max ~ 1: we are BHU1, matching the shallow hunch.",
        "Hotter, or an old supraverse, gives D_max >> 1 and the depth reverts to the",
        "branching ratio m -- deep. So 'are we BHU1-2' becomes a statement about the",
        "void temperature and the supraverse age, both genuinely unknown (Q15a/Q15c);",
        "our existence as n>=1 already forces T_bath above the sterile threshold.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "growth.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: depletion freeze-out vs runaway -- the Q15a regulation answer.
    tr, xr = g.simulate_runaway(x0=2.0, t_max=0.49)
    axL.plot(tr, xr, "C3", lw=2.0, label="static bath: runaway")
    for bath, col in [(20, "C0"), (80, "C1"), (300, "C2")]:
        td, xd, bd = g.simulate_depletion(x0=2.0, bath=bath, t_max=30.0)
        axL.plot(td, xd, col, lw=1.8, label=fr"deplete patch $b_0={bath}$")
        axL.axhline(2 + bath, color=col, ls=":", lw=0.9)
    axL.axhline(1.0, color="0.6", ls="--", lw=1.0)
    axL.text(8, 1.4, r"$M_{\rm crit}$", fontsize=8, color="0.4")
    axL.set_yscale("log")
    axL.set_xlabel(r"time  $t'=t/\tau_0$")
    axL.set_ylabel(r"mass  $M/M_{\rm crit}$")
    axL.set_title("Growth runs away, depletion freezes it\n"
                  r"$M_{\rm final}\simeq$ patch mass (regulates Q15a)", fontsize=11)
    axL.set_ylim(0.8, 1e3)
    axL.set_xlim(0, 30)
    axL.legend(fontsize=8, loc="center right")
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: depth cap vs void bath temperature -- the population payoff.
    T = np.logspace(10.5, 15.5, 400)
    for label, age, col in [("1 Hubble time", HUBBLE_TIME_S, "C0"),
                            (r"$10^{6}$ Hubble times", 1e6 * HUBBLE_TIME_S, "C1")]:
        D = np.array([max(g.max_depth(age, Ti), 0.1) for Ti in T])
        axR.plot(T, D, col, lw=2.0, label=f"age = {label}")
    axR.axhspan(0.1, 1.0, color="0.85", alpha=0.6)
    axR.text(2e11, 0.3, "sterile foam (no BHU) -- excluded by $n\\geq1$",
             fontsize=8, color="0.35")
    axR.axhspan(1.0, 2.0, color="C0", alpha=0.15)
    axR.text(2e11, 1.35, "BHU1-2 (forced shallow)", fontsize=8, color="C0")
    axR.text(2e11, 1e4, "deep: depth set by branching $m$", fontsize=8, color="0.35")
    axR.axhline(1.0, color="0.5", ls=":", lw=1.0)
    axR.axhline(2.0, color="0.5", ls=":", lw=1.0)
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"void bath temperature  $T_{\rm bath}$  [K]")
    axR.set_ylabel(r"max generation depth  $D_{\max}=T_{\rm supra}/\tau_{\rm gen}$")
    axR.set_title("The void temperature dials our depth\n"
                  r"$\tau_{\rm gen}\propto T_{\rm bath}^{-3}$", fontsize=11)
    axR.set_ylim(0.1, 1e9)
    axR.legend(fontsize=8.5, loc="lower right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "growth.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
