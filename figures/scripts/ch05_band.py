#!/usr/bin/env python3
r"""The viable band is narrow -> the branching ratio is small -> we are shallow.

Viable children need near-Hubble-mass progenitor holes, so mass conservation caps
the branching at m = epsilon f_clean (M_parent/M_vis). If the viable band is narrow
(M_vis ~ our mass, w = M_parent/M_vis ~ O(1)) then m ~ O(0.1) < 1: subcritical,
shallow, most likely BHU1. Only a wide band (w >> 1) gives a deep population.
Renders figures/pdf/viable_band.pdf, writes sims/output/viable_band.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import population as pop

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

EPS = 0.1            # fraction of parent mass reaching viable-size holes
F_CLEAN = 1.0        # big holes are fair-sample (clean-family) accretors


def main() -> None:
    lines = [
        "The viable band is narrow -> small branching -> shallow",
        "=" * 58,
        "  m = epsilon * f_clean * (M_parent/M_vis), capped at M_parent/M_vis.",
        f"  (epsilon={EPS}, f_clean={F_CLEAN}); w = M_parent/M_vis = size-band width.",
        "",
        "   w=M_p/M_vis     m       P(BHU1)  P(BHU2)  P(BHU1 or 2)",
    ]
    for w in (1.0, 3.0, 10.0, 1e2, 1e4, 1e6):
        m = pop.mass_budget_branching(w, 1.0, EPS, F_CLEAN)
        p1, p2 = pop.prob_bhu(m, 1), pop.prob_bhu(m, 2)
        lines.append(f"   {w:9.0e}    {m:7.2f}   {p1:6.3f}   {p2:6.3f}   "
                     f"{pop.shallow_probability(m):6.3f}")
    lines += [
        "",
        "READING: the old 'm ~ 1e18, so we are deep' estimate double-counted -- it",
        "treated every stellar hole (which makes a tiny ~1e31 kg, sub-viable child)",
        "as a universe. Viable children need ~Hubble-mass holes, and mass conserva-",
        "tion allows only ~M_parent/M_vis of those. If the viable band is narrow",
        "(M_vis ~ our mass, w ~ O(1) -- the edge we appear to sit at), then m ~ 0.1",
        "< 1: subcritical, shallow, most likely BHU1. The same big holes that meet",
        "the size cut are the fair-sample (clean-family) accretors, so the size and",
        "eta cuts select together -- the band is narrow on both axes at once. A wide",
        "band (M_vis << M_parent) would be needed for a deep population.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "viable_band.txt"), "w") as f:
        f.write(text + "\n")

    w = np.logspace(0.0, 6.0, 400)
    m = np.array([pop.mass_budget_branching(wi, 1.0, EPS, F_CLEAN) for wi in w])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # Panel L: branching ratio vs size-band width.
    axL.plot(w, m, "C0", lw=2.2)
    axL.axhline(1.0, color="0.5", ls=":", lw=1.1)
    axL.text(1.3, 1.25, "critical $m=1$", fontsize=8, color="0.4")
    w_crit = 1.0 / EPS / F_CLEAN
    axL.axvspan(1.0, w_crit, color="C0", alpha=0.13)
    axL.text(1.15, 30, "narrow band\n$m<1$: shallow\n(BHU1-2)", fontsize=8,
             color="C0")
    axL.text(200, 30, "wide band\n$m>1$: deep", fontsize=8, color="0.4")
    axL.plot([3.0], [pop.mass_budget_branching(3.0, 1.0, EPS, F_CLEAN)], "ko",
             ms=5)
    axL.annotate("us? (viable ~ our size)", xy=(3.0, 0.3), xytext=(4, 3),
                 fontsize=8, arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"viable size-band width  $w=M_{\rm parent}/M_{\rm vis}$")
    axL.set_ylabel(r"branching ratio  $m=\epsilon\,f_{\rm clean}\,w$")
    axL.set_title("Mass conservation caps the branching\n"
                  r"$m\leq M_{\rm parent}/M_{\rm vis}$ (not $\sim10^{18}$)",
                  fontsize=11)
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: resulting depth probabilities vs band width.
    for n, col in [(1, "C0"), (2, "C1"), (3, "C2")]:
        axR.plot(w, [pop.prob_bhu(pop.mass_budget_branching(wi, 1.0, EPS, F_CLEAN), n)
                     for wi in w], col, lw=2.0, label=fr"$P(\mathrm{{BHU}}_{n})$")
    axR.plot(w, [pop.shallow_probability(pop.mass_budget_branching(wi, 1.0, EPS, F_CLEAN))
                 for wi in w], "k--", lw=1.4, label=r"$P(\mathrm{BHU1\ or\ 2})$")
    axR.axvspan(1.0, w_crit, color="C0", alpha=0.10)
    axR.set_xscale("log")
    axR.set_xlabel(r"viable size-band width  $w=M_{\rm parent}/M_{\rm vis}$")
    axR.set_ylabel(r"$P(\,\cdot\mid n\geq1)$")
    axR.set_title("Narrow band $\\Rightarrow$ shallow; wide band $\\Rightarrow$ deep",
                  fontsize=11)
    axR.set_ylim(0, 1.02)
    axR.legend(fontsize=8.5, loc="center right")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "viable_band.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
