#!/usr/bin/env python3
r"""What a hole eats sets the child's eta; and are we likely BHU1 or BHU2?

Left: adiabatic bounce (D=1) => eta_child = eta of the accreted material. Our
photon-dominated eta ~ eta_cosmic places our progenitor at a fair-sample (b~1)
accretor, not a concentrated-baryon (star/galaxy) one -> two families of
universes. Right: given the floor n>=1 (we are not an OGU), P(BHU_n) ~ m^n, so
'BHU1 or 2' is likely only for subcritical branching m < 1 (or a young
supraverse). Renders figures/pdf/accretion_depth.pdf, writes output/accretion.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import accretion as ac
from cartasis_sims import population as pop

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    lines = [
        "What the hole eats sets the child's eta (adiabatic D=1 bounce)",
        "=" * 60,
        "  eta_child = b * eta_cosmic, b = eta_infall/eta_cosmic, capped at 1.",
        "",
        "  source                          b         eta_child    family",
    ]
    for name, b in ac.SOURCES.items():
        ec = ac.eta_child(b)
        fam = "baryon-rich" if ac.is_baryon_rich(b) else "clean"
        lines.append(f"  {name:30s} {b:7.1e}   {ec:.2e}   {fam}")
    lines += [
        "",
        "  Our eta ~ 6e-10 ~ eta_cosmic => b ~ 1: our progenitor accreted a FAIR",
        "  cosmic sample (a horizon-scale / radiation-fed hole), NOT concentrated",
        "  baryons. A star-fed hole (b~1e9) would have spawned a baryon-rich,",
        "  degenerate universe (eta ~ 1) -- a different family, not ours.",
        "",
        "Are we BHU1 or BHU2?  P(BHU_n | n>=1) ~ m^n (m = viable children/universe)",
        "-" * 60,
        "   m      P(BHU1)   P(BHU2)   P(n>=3)   P(BHU1 or 2)",
    ]
    for m in (0.1, 0.3, 0.5, 0.9, 1.0, 1.5, 3.0):
        p1, p2 = pop.prob_bhu(m, 1), pop.prob_bhu(m, 2)
        lines.append(f"  {m:4.1f}    {p1:6.3f}    {p2:6.3f}    {1-p1-p2:6.3f}    "
                     f"{pop.shallow_probability(m, n_max=2):6.3f}")
    lines += [
        "",
        "READING: shallow (BHU1-2) is likely ONLY if the branching ratio is",
        "subcritical, m < 1 -- on average <=1 viable child per universe. But the",
        "naive count is huge: a universe forms ~1e18+ black holes, ~all super-",
        "critical, ~half past the spin-survival threshold, so m ~ 1e18 unless",
        "viability is far more restrictive (e.g. only fair-sample / clean-fed holes",
        "spawn universes like ours, the rest making the baryon-rich family). With",
        "m >> 1 the population is dominated by the deepest generation the supraverse",
        "has had TIME to grow, and we are almost surely DEEP, not BHU1-2. Being",
        "shallow then requires either restrictive viability (m -> <1) or a YOUNG",
        "supraverse (small truncation D). This is the live tension; the DM+DE anchor",
        "(n>=1) is the only piece the data settle cleanly.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "accretion.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: eta_child vs accretion bias b -- the two families.
    b = np.logspace(-1.5, 9.5, 500)
    ec = np.array([ac.eta_child(bi) for bi in b])
    axL.plot(b, ec, "C0", lw=2.2)
    axL.axhline(ac.ETA_COSMIC, color="C3", lw=1.3)
    axL.text(2e-2, ac.ETA_COSMIC * 1.4, r"our $\eta\simeq6\times10^{-10}$",
             color="C3", fontsize=9)
    axL.axhspan(1e-2, 1.5, color="C1", alpha=0.13)
    axL.text(2e4, 8e-2, "baryon-rich / degenerate\nfamily ($\\eta\\to1$)",
             fontsize=8, color="C1")
    axL.axhspan(1e-12, 1e-6, color="C0", alpha=0.10)
    axL.text(2e-2, 4e-12, "clean / photon-dominated family", fontsize=8,
             color="C0")
    marks = [("fair sample\n(=us, $b{\\sim}1$)", 1.0), ("galaxy ISM", 3e6),
             ("a star", 1e9), ("void/IGM", 0.1)]
    for label, bi in marks:
        axL.plot([bi], [ac.eta_child(bi)], "ko", ms=4)
        axL.annotate(label, xy=(bi, ac.eta_child(bi)), fontsize=7.5,
                     xytext=(0, 8), textcoords="offset points", ha="center")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"accretion baryon-bias  $b=\eta_{\rm infall}/\eta_{\rm cosmic}$")
    axL.set_ylabel(r"child asymmetry  $\eta_{\rm child}=b\,\eta_{\rm parent}$")
    axL.set_title(r"What the hole eats sets the child's $\eta$"
                  "\n(adiabatic $D{=}1$: $\\eta_{\\rm child}=\\eta_{\\rm infall}$)",
                  fontsize=11)
    axL.set_ylim(1e-12, 3.0)
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: P(BHU_n) vs branching ratio m.
    ms = np.linspace(0.02, 3.0, 400)
    for n, col in [(1, "C0"), (2, "C1"), (3, "C2")]:
        axR.plot(ms, [pop.prob_bhu(m, n) for m in ms], col, lw=2.0,
                 label=fr"$P(\mathrm{{BHU}}_{n})$")
    axR.plot(ms, [pop.shallow_probability(m, n_max=2) for m in ms], "k--", lw=1.4,
             label=r"$P(\mathrm{BHU1\ or\ 2})$")
    axR.axvline(1.0, color="0.5", ls=":", lw=1.1)
    axR.text(1.03, 0.55, "critical\n$m=1$", fontsize=8, color="0.4")
    axR.axvspan(0.02, 1.0, color="C0", alpha=0.08)
    axR.text(0.1, 0.05, "subcritical:\nshallow, BHU1-2", fontsize=8, color="C0")
    axR.text(1.5, 0.05, "supercritical:\ndeep ($n\\gg1$)", fontsize=8, color="0.4")
    axR.set_xlabel(r"branching ratio  $m$ = mean viable children per universe")
    axR.set_ylabel(r"$P(\,\cdot\mid n\geq1)$")
    axR.set_title("Are we BHU1 or BHU2?\nonly if branching is subcritical ($m<1$)",
                  fontsize=11)
    axR.set_ylim(0, 1.02)
    axR.legend(fontsize=8.5, loc="upper right")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "accretion_depth.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
