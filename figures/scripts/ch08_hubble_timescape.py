#!/usr/bin/env python3
r"""The Hubble tension from inhomogeneity: the timescape clock-rate effect.

Left: the dressed/bare Hubble ratio R(f_v) = (4 f_v^2 + f_v + 4)/(2(2+f_v)) -- the
factor by which a wall observer (us) over-reads the volume-average Hubble. R(0)=1
(homogeneous, no effect); it crosses the observed 8.4% early-vs-late gap at
f_v~0.47 and reaches 28% at the observed void fraction f_v0~0.76. So inhomogeneity
at the OBSERVED void fraction has the right sign and more than enough magnitude.
Right: the tracker void fraction f_v(t)=t/(t+b) -- ~1e-4 at recombination (CMB
untouched), rising to 0.76 today: the effect is purely late-time, which is exactly
the early-vs-late structure of the tension. Renders figures/pdf/hubble_timescape.pdf,
writes sims/output/hubble_timescape.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import hubble_timescape as ht

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    tau = ht.observed_tension_fraction()
    f_match = ht.void_fraction_to_match_tension()
    offset = ht.wall_vs_average_offset()
    f_rec = ht.void_fraction_at_recombination()
    bare = ht.H0_DRESSED_BESTFIT / ht.dressed_over_bare(ht.FV0_BESTFIT)

    lines = [
        "Hubble tension from inhomogeneity (timescape clock-rate effect)",
        "=" * 62,
        f"Observed gap:  H0_local={ht.H0_LOCAL} vs H0_CMB={ht.H0_CMB} "
        f"-> tau=(local-CMB)/CMB = {tau*100:.1f}%",
        "",
        "Dressed/bare Hubble ratio  R(f_v) = (4 f_v^2 + f_v + 4)/(2(2+f_v)):",
        f"  R(0)      = {ht.dressed_over_bare(0.0):.3f}   (homogeneous: no effect)",
        f"  R({f_match:.2f})   = {ht.dressed_over_bare(f_match):.3f}   "
        f"(void fraction that reproduces the {tau*100:.1f}% gap exactly)",
        f"  R(0.762)  = {ht.dressed_over_bare(0.762):.3f}   "
        f"(observed f_v0 -> {offset*100:.0f}% wall-vs-average excess)",
        "",
        f"=> wall/local Hubble EXCEEDS volume-average/global by {offset*100:.0f}%,",
        f"   against an {tau*100:.1f}% tension: right sign, ample magnitude, and the",
        f"   void fraction needed ({f_match:.2f}) is BELOW the observed 0.76.",
        "",
        "Late-time only (the early-vs-late structure):",
        f"  f_v at recombination (z~1100) = {f_rec:.2e}  ->  R-1 = "
        f"{ht.dressed_over_bare(f_rec)-1:.1e}",
        "  i.e. the CMB sound-horizon H0 is untouched; the effect is purely local.",
        "",
        f"Timescape best fit (Wiltshire): dressed H0={ht.H0_DRESSED_BESTFIT}, "
        f"bare H0={bare:.1f}, gamma-bar0={ht.GAMMA0_BESTFIT}, f_v0={ht.FV0_BESTFIT}.",
        "",
        "SCT is natively inhomogeneous (void-dominated foam, parent-accretion dark",
        "sector, Ch6), so it inherits this resolution with NO extra parameter. This",
        "is a sign-and-magnitude result, not a substitute for the global data fit.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "hubble_timescape.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.8))

    # --- Panel 1: R(f_v) with the observed-tension band ---
    fv = np.linspace(0.0, 1.0, 400)
    R = np.array([ht.dressed_over_bare(x) for x in fv])
    axL.plot(fv, (R - 1.0) * 100.0, color="0.15", lw=2.2,
             label=r"timescape $R(f_v)-1$")
    axL.axhspan(0.0, tau * 100.0, color="0.55", alpha=0.18, lw=0,
                label=f"observed gap (0--{tau*100:.1f}%)")
    axL.axhline(tau * 100.0, color="0.45", ls="--", lw=1.0)
    # observed void fraction marker
    axL.axvline(ht.FV0_BESTFIT, color="0.30", ls=":", lw=1.0)
    axL.plot([ht.FV0_BESTFIT], [(ht.dressed_over_bare(ht.FV0_BESTFIT) - 1) * 100],
             "o", color="0.0", ms=7, zorder=6,
             label=f"observed $f_{{v0}}={ht.FV0_BESTFIT}$ ({offset*100:.0f}%)")
    axL.plot([f_match], [tau * 100.0], "s", color="0.45", ms=7, zorder=6,
             label=f"matches gap at $f_v={f_match:.2f}$")
    axL.set_xlabel(r"void volume fraction $f_v$")
    axL.set_ylabel(r"wall-vs-average Hubble excess $(R-1)$  [%]")
    axL.set_title("Inhomogeneity over-reads the local Hubble")
    axL.set_xlim(0, 1)
    axL.set_ylim(-2, 52)
    axL.legend(fontsize=7.5, loc="upper left")
    axL.grid(True, alpha=0.2)

    # --- Panel 2: void-fraction history, effect is late-time ---
    t = np.logspace(-5, 0.3, 400)        # t/t0
    fvt = np.array([ht.void_fraction(x) for x in t])
    axR.semilogx(t, fvt, color="0.15", lw=2.2, label=r"tracker $f_v(t)=t/(t+b)$")
    axR.axhline(ht.FV0_BESTFIT, color="0.45", ls="--", lw=1.0,
                label=f"today $f_{{v0}}={ht.FV0_BESTFIT}$")
    # recombination marker
    t_rec = (1.0 + 1100.0) ** (-1.5)
    axR.axvline(t_rec, color="0.30", ls=":", lw=1.0)
    axR.annotate("recombination\n$f_v\\sim10^{-4}$ (CMB untouched)",
                 xy=(t_rec, 0.02), xytext=(t_rec * 4, 0.30), fontsize=8,
                 color="0.2",
                 arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8))
    axR.plot([1.0], [ht.FV0_BESTFIT], "o", color="0.0", ms=7, zorder=6)
    axR.set_xlabel(r"time  $t/t_0$")
    axR.set_ylabel(r"void volume fraction $f_v$")
    axR.set_title("The clock effect is purely late-time")
    axR.set_ylim(0, 0.85)
    axR.legend(fontsize=8, loc="upper left")
    axR.grid(True, alpha=0.2, which="both")

    fig.suptitle("Hubble tension from inhomogeneity: the timescape clock-rate effect",
                 fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_pdf = os.path.join(PDF_DIR, "hubble_timescape.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
