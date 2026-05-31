#!/usr/bin/env python3
r"""OGU mass and spin distributions, from birth vs growth and birth vs viability.

Left: the OGU mass distribution. Births favour small masses (exp(-M/M0));
runaway growth (Mdot ~ M^g) transports mass upward, giving a stationary power law
n(M) ~ M^{-g} between the birth scale and a death cutoff -- many small, few large.
A low-mass viability cut (enough mass for astrophysics and black-hole formation)
selects observer-bearing universes, which pile up just above the edge: we sit
near the viability optimum, not in the rare massive tail.

Right: the OG spin distribution. A random seed most likely has LOW net vorticity
(Gaussian at 0), but below omega_min baryogenesis fails -- matter annihilates,
the universe is sterile and evaporates (no descendants). Viable, observer-bearing
seeds therefore sit just above threshold: low spin, low purity, just past hellish.
Our small observed asymmetry (eta ~ 1e-9..1e-8) places us right there -- testable
against any net rotation of our universe.

Writes figures/pdf/ogu_distributions.pdf.
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


def main() -> None:
    # ---- mass distribution ----
    M = np.logspace(-0.5, 3.2, 600)
    M0, M_vis, M_max = 1.0, 12.0, 800.0
    n_bondi = pop.ogu_mass_density(M, g=2.0, M0=M0, M_max=M_max)
    n_edd = pop.ogu_mass_density(M, g=1.0, M0=M0, M_max=M_max)
    n_obs = pop.observer_mass_density(M, g=2.0, M0=M0, M_max=M_max, M_vis=M_vis)
    n_bondi /= n_bondi.max()
    n_edd /= n_edd.max()
    n_obs /= n_obs.max() if n_obs.max() > 0 else 1.0
    M_us = M[np.argmax(n_obs)]

    # ---- spin distribution ----
    w = np.linspace(0, 5, 1500)
    sigma, wmin = 1.0, 0.6
    p_birth = pop.spin_birth_pdf(w, sigma)
    p_obs = pop.spin_observer_pdf(w, sigma, omega_min=wmin, prod_power=1.0)
    p_via = np.where(w > wmin, p_birth, 0.0)       # viability-only (no productivity)
    p_birth /= p_birth.max()
    p_obs /= p_obs.max() if p_obs.max() > 0 else 1.0
    p_via /= p_via.max() if p_via.max() > 0 else 1.0
    w_peak = w[np.argmax(p_obs)]
    sterile = pop.sterile_fraction(sigma, wmin)

    lines = [
        "OGU mass and spin distributions (birth vs growth / viability)",
        "=" * 60,
        "MASS: births ~exp(-M/M0); growth Mdot~M^g => n(M)~M^{-g} (power law).",
        f"  Bondi g=2 and Eddington g=1 shown; viability edge M_vis={M_vis}.",
        f"  Observer-bearing universes peak at M ~ {M_us:.1f} (just above M_vis):",
        "  we sit near the viability optimum, not the rare massive tail.",
        "",
        "SPIN: P_birth Gaussian at 0 (low spin likeliest); baryogenesis needs",
        f"  |omega|>omega_min={wmin} (eta = C|omega|/T > eta_min), else sterile.",
        f"  Sterile (sub-threshold, evaporating) fraction of seeds: {sterile:.2f}.",
        f"  Viable/observer spin peaks just above threshold at omega ~ {w_peak:.2f}:",
        "  viable universes are LOW-purity, just past hellish.  Our small eta",
        "  (~1e-9..1e-8) places us at that low-spin edge -- compare to any net",
        "  rotation of our universe (galaxy-spin handedness, CMB axis hints).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "ogu_distributions.txt"), "w") as f:
        f.write(text + "\n")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.7))

    ax1.loglog(M, n_bondi, "C0", lw=1.8, label="all OGUs, Bondi $n\\propto M^{-2}$")
    ax1.loglog(M, n_edd, "C0", lw=1.2, ls="--",
               label="all OGUs, Eddington $n\\propto M^{-1}$")
    ax1.fill_between(M, 1e-6, n_obs, where=(n_obs > 0), color="C2", alpha=0.25,
                     label="observer-bearing (viable)")
    ax1.axvline(M_vis, color="C3", lw=1.0, ls=":")
    ax1.annotate("viability edge\n(astrophysics, BHs)", xy=(M_vis, 3e-3),
                 xytext=(M_vis * 1.4, 2e-2), fontsize=7.5, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    ax1.plot([M_us], [1.0], "kv", ms=9, label="$\\sim$ our universe")
    ax1.set_xlabel("OGU mass $M$  (birth-scale units)")
    ax1.set_ylabel("relative number $n(M)$")
    ax1.set_ylim(1e-5, 2)
    ax1.set_title("Mass: birth (small) vs growth (large)")
    ax1.legend(fontsize=7.5, loc="lower left")
    ax1.grid(True, which="both", alpha=0.15)

    ax2.plot(w, p_birth, "0.5", lw=1.5, label="birth $P(\\omega)$ (low spin)")
    ax2.plot(w, p_via, "C0", lw=1.3, ls="--", label="viable only")
    ax2.fill_between(w, 0, p_obs, color="C2", alpha=0.30,
                     label="observer-bearing ($\\times$ productivity)")
    ax2.axvspan(0, wmin, color="C3", alpha=0.12)
    ax2.annotate("sterile / hellish\n(annihilates, evaporates)",
                 xy=(wmin * 0.5, 0.5), xytext=(wmin * 0.15, 0.62),
                 fontsize=7.5, color="C3")
    ax2.axvline(wmin, color="C3", lw=1.0, ls=":")
    ax2.plot([w_peak], [1.0], "kv", ms=9, label="$\\sim$ our universe (small $\\eta$)")
    ax2.set_xlabel("OG seed net vorticity $|\\omega|$  ($\\sigma$ units)")
    ax2.set_ylabel("relative probability")
    ax2.set_title("Spin: low likeliest, but needs $\\omega>\\omega_{\\min}$ to live")
    ax2.legend(fontsize=7.5, loc="upper right")
    ax2.grid(True, alpha=0.15)

    fig.suptitle("OGU populations: where birth, growth, and viability put us",
                 fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(PDF_DIR, "ogu_distributions.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
