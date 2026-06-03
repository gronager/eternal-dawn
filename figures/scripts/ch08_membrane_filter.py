#!/usr/bin/env python3
r"""The membrane filter and the dark-to-baryon ratio f ~ 1/6 (Ch.8, Tier 3).

If dark matter is the parent matter the bounce membrane rejects, the dark-to-baryon
ratio is the membrane PASS-FRACTION f = omega_b/(omega_b+omega_c) ~ 0.157 ~ 1/6.
Modelling the membrane as a torsion barrier, the pass-fraction is a penetration /
occupation factor in x = Delta/T (barrier-to-thermal ratio): f = exp(-x) (Boltzmann)
or f = 1/(exp(x)+1) (Fermi). The membrane is, by definition, where the torsion energy
equals the matter energy, so x ~ O(1) -- giving f in the 0.08-0.37 decade, with the
observed 1/6 at x ~ 1.8. This NARROWS f from a free parameter to an O(1) barrier the
bounce supplies; it does not yet pin 0.157 to three digits (that needs the full bounce
profile). Renders figures/pdf/membrane_filter.pdf, writes sims/output/membrane_filter.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import membrane_filter as mf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    f_obs = mf.observed_f()
    xB = mf.barrier_ratio_for_f(f_obs, "boltzmann")
    xF = mf.barrier_ratio_for_f(f_obs, "fermi")
    bandB = mf.natural_band("boltzmann")
    bandF = mf.natural_band("fermi")

    lines = [
        "The membrane filter and the dark-to-baryon ratio f ~ 1/6",
        "=" * 58,
        f"  observed pass-fraction f = omega_b/(omega_b+omega_c) = {f_obs:.4f} "
        f"(~1/6; omega_c/omega_b ~ 5.4)",
        f"  bounce thermal scale T ~ {mf.BOUNCE_T_GEV:.1e} GeV",
        "",
        "  model: membrane = torsion barrier; pass-fraction f(x), x = Delta/T",
        f"    f = 1/6  <=>  x = ln6 = {np.log(6):.3f} (Boltzmann),  "
        f"x = ln5 = {np.log(5):.3f} (Fermi)",
        f"    observed f sits at  x = {xB:.2f} (Boltzmann) / {xF:.2f} (Fermi)",
        "",
        "  anchor: at the membrane torsion energy = matter energy (bounce balance),",
        "          so x = Delta/T ~ O(1) -- not 1e-9, not 1e9.",
        f"    an O(1) barrier x in [{mf.X_NATURAL_LO}, {mf.X_NATURAL_HI}] gives "
        f"f in [{bandB[0]:.3f}, {bandB[1]:.3f}] (Boltzmann)",
        f"    / [{bandF[0]:.3f}, {bandF[1]:.3f}] (Fermi) -- the observed 1/6 is inside.",
        "",
        "READING: f ~ 1/6 is NATURAL -- the generic output of an O(1) torsion barrier,",
        "not a tuned number. What is owed is the precise x ~ 1.8 (vs 1.0 or 2.5), which",
        "needs the full inhomogeneous bounce profile (Tier-3). So f goes from a free fit",
        "in [0,1] to an O(1) ratio the bounce supplies -- narrowed, not yet pinned.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "membrane_filter.txt"), "w") as fobj:
        fobj.write(text + "\n")

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    x = np.linspace(0.0, 4.0, 400)
    ax.plot(x, np.exp(-x), "C0", lw=2.2, label=r"Boltzmann  $f=e^{-x}$")
    ax.plot(x, 1.0 / (np.exp(x) + 1.0), "C1", lw=2.2,
            label=r"Fermi  $f=1/(e^{x}+1)$")

    # the O(1)-barrier band (bounce balance)
    ax.axvspan(mf.X_NATURAL_LO, mf.X_NATURAL_HI, color="0.85", alpha=0.7)
    ax.text(1.75, 0.52, "O(1) torsion barrier\n(bounce balance: "
            r"$\Delta\!\sim\!T$)", ha="center", fontsize=9, color="0.3")

    # the observed f and where it lands
    ax.axhline(f_obs, color="C3", ls="--", lw=1.5)
    ax.text(3.4, f_obs + 0.012, rf"observed $f={f_obs:.3f}\approx1/6$",
            color="C3", fontsize=9.5, ha="right")
    for xv, c in ((xB, "C0"), (xF, "C1")):
        ax.plot(xv, f_obs, "o", color=c, ms=8)
        ax.annotate(rf"$x={xv:.2f}$", (xv, f_obs), textcoords="offset points",
                    xytext=(2, -16), fontsize=8.5, color=c)

    ax.set_xlim(0, 4)
    ax.set_ylim(0, 0.6)
    ax.set_xlabel(r"barrier-to-thermal ratio  $x=\Delta/T$")
    ax.set_ylabel(r"membrane pass-fraction  $f=\omega_b/(\omega_b+\omega_c)$")
    ax.set_title("The dark-to-baryon ratio as a torsion-barrier pass-fraction\n"
                 r"$f\sim1/6$ is the natural output of an $O(1)$ barrier, not a tuned number",
                 fontsize=11.5)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "membrane_filter.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "membrane_filter.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'membrane_filter.pdf')}")


if __name__ == "__main__":
    main()
