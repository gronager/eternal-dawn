#!/usr/bin/env python3
r"""The dark sector as look-back into the parent: gravity felt before light seen.

The membrane is our PAST horizon. Parent-injected matter arrives in our deep past, so
(left) the particle horizon grows every year -- we see further back with time, toward
a de Sitter asymptote -- while (right) we feel the parent's gravity (dark matter) NOW,
ahead of ever seeing its light: the dark sector is a rate, not a relic.
Renders figures/pdf/lookback.pdf, writes sims/output/lookback.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import lookback as lb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    h_now = lb.horizon_gly(1.0)
    h_inf = lb.asymptotic_horizon_gly()
    gain = lb.lookback_gain(1.0)
    frac = lb.fraction_seen_now()
    t0 = lb.cosmic_age(1.0) / lb.GYR_S

    lines = [
        "The dark sector as look-back into the parent",
        "=" * 58,
        f"  cosmic age now                = {t0:.2f} Gyr",
        f"  comoving particle horizon now = {h_now:.1f} Gly  (how far back we see)",
        f"  asymptotic horizon (a->inf)   = {h_inf:.1f} Gly  (the most we will EVER see)",
        f"  fraction of that seen today   = {frac:.2f}",
        f"  look-back gain after +1 Gyr   = +{gain:.2f} Gly  (we see FURTHER back)",
        "",
        "FULL DYNAMIC (all three, consistent):",
        "  1. We see further back every year. The past (particle) horizon recedes",
        "     monotonically toward ~64 Gly, so parent matter injected at the membrane",
        "     becomes visible only as the horizon expands out to it. Look in Hubble",
        "     today: ~14 Gyr back. Look again in +1 Gyr: ~1 Gly further back -- yes.",
        "  2. We feel its gravity before we see it. The parent's projected mass is",
        "     already in our metric (it curves our space NOW = dark matter), while its",
        "     light is still climbing in from the past horizon. Gravity leads light.",
        "  3. So the dark sector is a RATE: dark matter ~ the parent's accumulated",
        "     projected mass, dark energy ~ the time-derivative of its growth. Both",
        "     evolve with the parent's meal -- the handle on sigma_8/S_8 and DESI w0-wa.",
        "",
        "READING: 'add more matter to a BHU and it arrives further in our past' is",
        "exactly right -- the membrane is our Big-Bang surface, so new parent matter",
        "lands at the receding past horizon, felt gravitationally before it is ever",
        "seen. Dark matter is the parent's pull arriving ahead of its light.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "lookback.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel L: particle horizon vs cosmic age, toward the de Sitter asymptote.
    a = np.logspace(-2, 1.2, 200)
    ages = np.array([lb.cosmic_age(ai) / lb.GYR_S for ai in a])
    hor = np.array([lb.horizon_gly(ai) for ai in a])
    axL.plot(ages, hor, "C0", lw=2.4)
    axL.axhline(h_inf, color="0.5", ls="--", lw=1.2)
    axL.text(2, h_inf + 1.0, f"de Sitter asymptote ~{h_inf:.0f} Gly "
             "(all we will ever see)", fontsize=8, color="0.4")
    axL.plot([t0], [h_now], "ko", ms=7)
    axL.annotate(f"now: {t0:.1f} Gyr, see {h_now:.0f} Gly back\n"
                 f"+1 Gyr $\\Rightarrow$ +{gain:.1f} Gly further",
                 xy=(t0, h_now), xytext=(15, 25), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.set_xlabel("cosmic age  [Gyr]")
    axL.set_ylabel("comoving particle horizon  [Gly]")
    axL.set_title("We see further back every year\n"
                  "(the past horizon recedes toward the de Sitter limit)", fontsize=11)
    axL.set_xlim(0, 40)
    axL.set_ylim(0, 70)
    axL.grid(True, alpha=0.2)

    # Panel R: gravity (felt now) leads light (seen later).
    r = np.linspace(0, 1, 200)
    seen = np.where(r <= frac, 1.0, 0.0)                 # visible (light arrived)
    axR.fill_between(r, 0, np.where(r <= frac, 1.0, 0.0), step="pre", color="C0",
                     alpha=0.25, label="seen (light arrived)")
    axR.fill_between(r, 0, 1.0, color="C3", alpha=0.12,
                     label="felt (gravity already in our metric)")
    axR.axvline(frac, color="C0", lw=1.6)
    axR.text(frac - 0.02, 0.5, "light\nhorizon", fontsize=8.5, color="C0", ha="right")
    axR.text(0.985, 0.5, "membrane\n(parent matter:\nfelt, not yet seen)",
             fontsize=8.5, color="C3", ha="right")
    axR.annotate("dark matter = the parent's gravity\narriving AHEAD of its light",
                 xy=(0.85, 0.78), xytext=(0.30, 0.84), fontsize=9, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3"))
    axR.set_xlabel("fraction of the way to the membrane (our past horizon)")
    axR.set_yticks([])
    axR.set_title("We feel the parent's mass before we see it\n"
                  "$\\Rightarrow$ dark matter is a rate, not a relic", fontsize=11)
    axR.set_xlim(0, 1)
    axR.set_ylim(0, 1)
    axR.legend(fontsize=8.5, loc="lower left")

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "lookback.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
