#!/usr/bin/env python3
r"""OGU genesis: the spin window (the pump) and the entrained mass (horizon mass).

Left: M_OGU ~ c^3 t/(2G) -- the rotating bounce entrains the void within its
causal reach until the patch fills and collapses to a horizon-scale hole, so the
OGU mass measures the supraverse age (our mass at ~1 Hubble time, ~1e65 kg at
~1e12 Hubble times -> ~1e11 siblings). Right: the spin window -- a seed must clear
the survival floor f_min (else sterile) but is choked near extremal, so the vortex
pump is most efficient at an intermediate spin sweet spot.
Renders figures/pdf/genesis.pdf, writes sims/output/genesis.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import genesis as g

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    fmin = g.spin_floor()
    fsweet = 1.0 / np.sqrt(2.0)

    lines = [
        "OGU genesis: vortex pump (spin window) and entrained mass (horizon mass)",
        "=" * 60,
        "  M_OGU ~ c^3 t/(2G) for a condensed void (patch fills -> horizon hole):",
        "    supraverse age      M_OGU          BHU1 siblings (eps=0.1)",
    ]
    for mult in (1.0, 1e3, 1e6, 1e9, 1e12):
        t = mult * g.HUBBLE_TIME_S
        M = g.ogu_mass(t)
        lines.append(f"    {mult:.0e} Hubble t   {M:.2e} kg     {g.n_siblings(M):.2e}")
    lines += [
        "",
        f"  -> our mass (~9e52 kg) = horizon mass at ~1 Hubble time;",
        f"     1e65 kg OGU <-> age ~{g.age_for_mass(1e65)/g.HUBBLE_TIME_S:.0e} "
        f"Hubble times <-> ~{g.n_siblings(1e65):.0e} sibling universes.",
        "",
        "  Spin window (eta_min=2e-11, C=10): the pump",
        f"    survival floor f_min = {fmin:.3f} (below: sterile, pumps nothing)",
        f"    sweet spot     f ~ {fsweet:.3f} (frame-drag vs centrifugal choking)",
        f"    near extremal  f -> 1: choked (matter circularizes, no radial infall)",
        "",
        "READING: the rotating bounce drags frames and entrains the surrounding",
        "void -- a gravitational vortex feeding the new universe. It fills its",
        "causal patch and collapses to a horizon-scale hole, so M_OGU ~ c^3 t/(2G):",
        "the OGU mass measures the supraverse age, and is DECOUPLED from our being",
        "BHU1 (that is fixed by descendants pinned at M_vis). The void density only",
        "matters if the void is sub-critical (then M_OGU is smaller). The spin must",
        "clear the survival floor but avoid centrifugal choking, so genesis favours",
        "an intermediate sweet spot -- the productivity factor behind the observed",
        "low-spin, low-purity viable universes (Ch.5 spin distribution).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "genesis.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: M_OGU vs supraverse age (horizon mass), with sibling-count axis.
    mult = np.logspace(-0.5, 13.0, 400)
    t = mult * g.HUBBLE_TIME_S
    M = np.array([g.ogu_mass(ti) for ti in t])
    axL.plot(mult, M, "C0", lw=2.2, label=r"condensed void: $M\sim c^3t/2G$")
    # a sub-critical void (density-limited, smaller)
    Msub = np.array([g.ogu_mass(ti, 1e-4 * g.void_critical_density(ti)) for ti in t])
    axL.plot(mult, Msub, "C1", lw=1.6, ls="--",
             label=r"sub-critical void ($10^{-4}\rho_{\rm crit}$)")
    for lbl, mm in [("us / our mass\n(1 Hubble time)", 1.0),
                    (r"$10^{65}$ kg: $\sim10^{11}$ siblings", g.age_for_mass(1e65) / g.HUBBLE_TIME_S)]:
        axL.plot([mm], [g.ogu_mass(mm * g.HUBBLE_TIME_S)], "ko", ms=5)
        axL.annotate(lbl, xy=(mm, g.ogu_mass(mm * g.HUBBLE_TIME_S)),
                     xytext=(0, -36), textcoords="offset points", fontsize=8,
                     ha="center", arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"supraverse age  $T_{\rm supra}$  [Hubble times]")
    axL.set_ylabel(r"OGU mass  $M_{\rm OGU}$  [kg]")
    axL.set_title("The entrained mass is the horizon mass\n"
                  r"$M_{\rm OGU}\sim c^3 T_{\rm supra}/2G$ (decoupled from our depth)",
                  fontsize=11)
    axL.legend(fontsize=8.5, loc="upper left")
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: the spin window / vortex sweet spot.
    f = np.linspace(0.0, 1.0, 500)
    eff = np.array([g.entrainment_efficiency(fi) for fi in f])
    axR.plot(f, eff / eff.max(), "C0", lw=2.2, label="entrainment efficiency")
    axR.axvspan(0.0, fmin, color="0.85", alpha=0.7)
    axR.text(fmin / 2, 0.5, "sterile\n($\\eta<\\eta_{\\min}$)", fontsize=8,
             color="0.35", ha="center", rotation=90)
    axR.axvline(fsweet, color="C3", ls=":", lw=1.2)
    axR.text(fsweet + 0.02, 0.92, "sweet spot", fontsize=8, color="C3")
    axR.text(0.9, 0.3, "choked\n(extremal)", fontsize=8, color="0.4", ha="center")
    # eta(f) overlaid on a second axis
    ax2 = axR.twinx()
    ax2.plot(f, [g.eta_from_spin(fi) for fi in f], "C2", lw=1.4, ls="--")
    ax2.set_ylabel(r"seed asymmetry  $\eta=Cf\,\Omega/T$", color="C2", fontsize=9)
    ax2.tick_params(axis="y", labelcolor="C2")
    axR.set_xlabel(r"spin fraction  $f=\omega/\Omega_{\rm bounce}$")
    axR.set_ylabel("vortex pump efficiency (normalized)")
    axR.set_title("The pump has a spin window\n"
                  "survival floor below, centrifugal choking above", fontsize=11)
    axR.set_xlim(0, 1)
    axR.set_ylim(0, 1.05)
    axR.legend(fontsize=8.5, loc="upper right")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "genesis.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
