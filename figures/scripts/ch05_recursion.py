#!/usr/bin/env python3
r"""The OGU as a finite-lived recursive timeline: saturate, evaporate, recurse.

Left: temperatures -- an OGU grows until its Hawking temperature falls to the void's
de Sitter floor T_dS, where it SATURATES (thermal equilibrium with its own horizon) at
M_eq ~ a Hubble mass. Right: the timescale ladder -- grow (~1 Hubble time) << interior
recursion tick (~1e108 s) << membrane lifetime (~1e142 s, ~1e34 interior cycles) <<
neighbour distance (e^{I/4} horizons, never reached). Every node is mortal; the tree
is eternal. Renders figures/pdf/recursion.pdf, writes sims/output/recursion.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import recursion as rc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    M_eq = rc.equilibrium_mass()
    T_dS = rc.de_sitter_temperature()
    t_life = rc.membrane_lifetime()
    n_cyc = rc.recursion_cycles_per_lifetime()

    lines = [
        "The OGU as a finite-lived recursive timeline",
        "=" * 58,
        f"  saturation mass M_eq (T_BH = T_dS) = {M_eq:.2e} kg  (~a Hubble mass, ~ours)",
        f"  void floor temperature T_dS        = {T_dS:.2e} K",
        f"  membrane lifetime (t_evap at M_eq) = {t_life:.2e} s  (~1e142 s)",
        f"  interior recursion tick            = {rc.T_INTERIOR_TICK:.0e} s  (~1e100 yr)",
        f"  recursion cycles per OGU lifetime  = {n_cyc:.2e}  (~1e34)",
        "",
        "TIMESCALE LADDER (each step astronomically longer than the last):",
        "  grow to saturation  ~ 1 Hubble time   ~ 1e18 s",
        "  interior recursion tick               ~ 1e108 s",
        "  membrane lifetime (then evaporates)   ~ 1e142 s",
        "  distance to nearest neighbour         ~ e^(I/4) horizons (never reached)",
        "",
        "READING: the membrane does NOT last forever -- it evaporates after ~1e142 s,",
        "long before it could ever touch a neighbour (which is e^(1e84) horizons away).",
        "But in that lifetime its interior runs ~1e34 full recursive cycles, each",
        "seeding its own sub-OGUs. So 'recursive timelines forever, or pulsating in the",
        "flat void?' is BOTH: every membrane is MORTAL (no branch nests literally",
        "forever), yet new OGUs nucleate without end -- in the flat void AND inside every",
        "heat-dead interior -- so the tree is eternal and ever-expanding though every",
        "node dies. A fractal of finite, pulsating timelines. The void is the one truly",
        "eternal, flat substrate; everything in it is a mortal recursive bubble.",
        "When an OGU evaporates, its inside (already a heat-dead de Sitter void hosting",
        "its own scatter of sub-OGUs) simply rejoins the flat void it came from -- mass",
        "and information returned, the cycle closed.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "recursion.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: T_BH(M) crossing the de Sitter floor -> saturation.
    M = np.logspace(48, 54, 400)
    T_BH = np.array([rc.hawking_temperature(m) for m in M])
    axL.plot(M, T_BH, "C0", lw=2.2, label=r"$T_{\rm BH}\propto 1/M$ (Hawking)")
    axL.axhline(T_dS, color="C3", lw=1.6, ls="--",
                label=r"$T_{\rm dS}$ (void floor)")
    axL.axvline(M_eq, color="0.5", ls=":", lw=1.1)
    axL.plot([M_eq], [T_dS], "ko", ms=7)
    axL.annotate("saturation\n$M_{\\rm eq}\\sim$ Hubble mass\n$T_{\\rm BH}=T_{\\rm dS}$:"
                 " growth stops", xy=(M_eq, T_dS), xytext=(1.5e48, 3e-29),
                 fontsize=8.5, arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.axvspan(M_eq, 1e54, color="C3", alpha=0.08)
    axL.text(7e52, 2e-28, "equilibrium\nwith the void\n(evaporates)", fontsize=8,
             color="C3")
    axL.text(2e48, 5e-31, "grows\n(absorbs net)", fontsize=8, color="C0")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"OGU mass  $M$  [kg]")
    axL.set_ylabel(r"temperature  [K]")
    axL.set_title("An OGU grows until it equals its void's temperature\n"
                  "then saturates at a Hubble mass (dark-energy capped)", fontsize=11)
    axL.legend(fontsize=8.5, loc="upper right")
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: the timescale ladder.
    stages = [
        ("grow to\nsaturation", 1e18, "C0"),
        ("interior\nrecursion tick", 1e108, "C2"),
        ("membrane\nlifetime\n(evaporates)", 1e142, "C3"),
        ("neighbour\ndistance\n(never reached)", 1e200, "0.5"),
    ]
    ypos = np.arange(len(stages))[::-1]
    for (lbl, t, col), y in zip(stages, ypos):
        axR.barh(y, np.log10(t), color=col, alpha=0.55, height=0.55)
        axR.text(np.log10(t) + 2, y, lbl, va="center", fontsize=8.5, color="0.15")
        axR.text(2, y, f"$10^{{{int(np.log10(t))}}}$ s", va="center", fontsize=8,
                 color="white", fontweight="bold")
    axR.set_yticks([])
    axR.set_xlim(0, 230)
    axR.set_xlabel(r"$\log_{10}$ (timescale / s)")
    axR.set_title("Every node is mortal, the tree is eternal\n"
                  r"$\sim10^{34}$ interior cycles per OGU, then it returns to the void",
                  fontsize=11)
    axR.grid(True, axis="x", alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "recursion.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
