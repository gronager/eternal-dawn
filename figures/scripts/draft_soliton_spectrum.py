#!/usr/bin/env python3
r"""Gate calculation for Parts II/III: the gravity-torsion soliton spectrum and the
electroweak S parameter (first approximation).

A: the self-generated four-fermion well M(r)=m+S(r) with the computed Dirac energy
   levels (the soliton's internal/relative ladder) and the ground-state wavefunction.
B: m^2 = E_n^2 vs level n -- the relativistic confining tower is approximately Regge
   (E^2 linear in n); a linear scalar well is cleanest, the anharmonic bounce well
   bends up.
C: the Peskin-Takeuchi S parameter vs M_V/f_pi. A QCD-like sector sits at S~0.25 (the
   technicolor graveyard, above the LEP bound); S<0.1 needs M_V/f_pi > 13 ('walking',
   anharmonic) -- NOT delivered here, the make-or-break Part III still owes.

Renders figures/pdf/soliton_spectrum.pdf, writes sims/output/soliton_spectrum.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import soliton as so

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    depth, width = 6.0, 1.0
    levels = {k: so.energy_levels(n_levels=6, kind=k, depth=depth, width=width,
                                  E_max=26.0, n_scan=130)
              for k in ("linear", "harmonic", "bounce")}
    lv = levels["bounce"]

    S_qcd = so.s_parameter(so.QCD_FPI_OVER_MV, so.QCD_MA_OVER_MV)
    need = so.fpi_over_MV_for_S(0.1, so.QCD_MA_OVER_MV)

    lines = [
        "Gravity-torsion soliton spectrum + electroweak S (first approximation)",
        "=" * 70,
        "  the 'potential' is the fermion's own four-fermion (Hehl-Datta) self-",
        "  interaction -- a nonlinear, self-consistent Dirac problem; here solved as a",
        "  first pass with the well held fixed (the spectrum SHAPE, in dimensionless",
        "  units; the absolute scale is the open density question).",
        "",
        "  energy levels (bounce well):  E_n = " + np.array2string(np.round(lv, 2)),
    ]
    for k in ("linear", "harmonic", "bounce"):
        sl, ic, rms = so.regge_slope(levels[k])
        lines.append(f"    {k:9s}: E^2 vs n Regge fit slope={sl:5.1f}, frac.rms={rms:.3f}"
                     + ("  (cleanest Regge)" if k == "linear" else ""))
    lines += [
        "",
        f"  electroweak S parameter (Peskin-Takeuchi, Weinberg sum rules):",
        f"    QCD-like (f_pi/M_V={so.QCD_FPI_OVER_MV:.2f}, M_A/M_V={so.QCD_MA_OVER_MV:.2f})"
        f"  ->  S = {S_qcd:.2f}   (> {so.S_LEP_BOUND}: the graveyard)",
        f"    S < 0.1 needs f_pi/M_V < {need:.3f}, i.e. M_V/f_pi > {1/need:.1f}"
        " (vs QCD's ~8): 'walking'.",
        "",
        "VERDICT (honest): the spectrum is approximately Regge (good, matches hadrons).",
        "But S < 0.1 is NOT established -- the naive value is ~0.25 (dead), and escaping",
        "needs walking/anharmonic dynamics the bounce well MIGHT supply but which we have",
        "not computed (it requires the composite vector/axial spectrum, not the single-",
        "particle tower). This is Part III's make-or-break, still owed.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "soliton_spectrum.txt"), "w") as f:
        f.write(text + "\n")

    fig = plt.figure(figsize=(15.0, 4.8))
    axA = fig.add_subplot(1, 3, 1)
    axB = fig.add_subplot(1, 3, 2)
    axC = fig.add_subplot(1, 3, 3)

    # ---- A: the self-generated well + levels + ground state ----
    r = np.linspace(0, 2.2, 400)
    M = so.effective_mass_profile(r, m=1.0, kind="bounce", depth=depth, width=width)
    axA.plot(r, M, "C3", lw=2.2, label=r"$M(r)=m+S(r)$ (four-fermion well)")
    axA.fill_between(r, 1.0, M, color="C3", alpha=0.07)
    for i, E in enumerate(lv[:5]):
        rt = so._turning_point(E, 1.0, "bounce", depth, width)
        axA.hlines(E, 0, rt, color="0.4", lw=1.0)
        axA.text(rt + 0.03, E, f"$E_{i}$", fontsize=8, color="0.4", va="center")
    rg, Gg, Fg = so.wavefunction(lv[0], kind="bounce", depth=depth, width=width)
    Gg = Gg / np.max(np.abs(Gg)) * 2.0 + lv[0]
    axA.plot(rg, Gg, "C0", lw=1.2, alpha=0.8, label="ground-state $G(r)$")
    axA.set_xlim(0, 2.2); axA.set_ylim(0, 18)
    axA.set_xlabel(r"radius $r$ (soliton units)")
    axA.set_ylabel(r"energy / effective mass")
    axA.set_title("The self-generated well and its levels", fontsize=10.5)
    axA.legend(fontsize=8, loc="upper left")
    axA.grid(True, alpha=0.2)

    # ---- B: Regge m^2 vs n ----
    for k, c in [("linear", "C0"), ("harmonic", "C1"), ("bounce", "C3")]:
        L = levels[k]
        n = np.arange(len(L))
        axB.plot(n, L**2, "o-", color=c, lw=1.4, ms=5, label=k)
    sl, ic, _ = so.regge_slope(levels["linear"])
    nn = np.arange(len(levels["linear"]))
    axB.plot(nn, sl * nn + ic, "k--", lw=0.9, alpha=0.6, label="linear Regge fit")
    axB.set_xlabel("level index $n$")
    axB.set_ylabel(r"$m_n^2 = E_n^2$")
    axB.set_title("Relativistic tower is ~Regge ($m^2\\propto n$)", fontsize=10.5)
    axB.legend(fontsize=8.5, loc="upper left")
    axB.grid(True, alpha=0.2)

    # ---- C: the S parameter wall ----
    MV_over_fpi = np.linspace(5, 20, 200)
    S = np.array([so.s_parameter(1.0 / x, so.QCD_MA_OVER_MV) for x in MV_over_fpi])
    axC.plot(MV_over_fpi, S, "C2", lw=2.2)
    axC.axhline(so.S_LEP_BOUND, color="0.4", ls="--", lw=1.2)
    axC.text(5.3, so.S_LEP_BOUND + 0.01, "LEP bound $S\\approx0.1$", fontsize=8, color="0.35")
    axC.axhspan(so.S_LEP_BOUND, 0.4, color="0.85", alpha=0.7)
    axC.plot(1 / so.QCD_FPI_OVER_MV, S_qcd, "ro", ms=8)
    axC.annotate("QCD-like\n(graveyard)", (1 / so.QCD_FPI_OVER_MV, S_qcd),
                 textcoords="offset points", xytext=(6, 0), fontsize=8, color="r")
    axC.plot(1 / need, 0.1, "C0o", ms=8)
    axC.annotate("walking\n(needs $M_V/f_\\pi>13$)", (1 / need, 0.1),
                 textcoords="offset points", xytext=(6, 6), fontsize=8, color="C0")
    axC.set_xlim(5, 20); axC.set_ylim(0, 0.4)
    axC.set_xlabel(r"$M_V/f_\pi$  (resonance heaviness / 'walking')")
    axC.set_ylabel("Peskin-Takeuchi $S$")
    axC.set_title("S: the make-or-break (not yet delivered)", fontsize=10.5)
    axC.grid(True, alpha=0.2)

    fig.suptitle("Gate calculation: the gravity-torsion soliton tower (Regge), and the "
                 "electroweak $S$ it must survive", fontsize=12.5, y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "soliton_spectrum.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "soliton_spectrum.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'soliton_spectrum.pdf')}")


if __name__ == "__main__":
    main()
