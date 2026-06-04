#!/usr/bin/env python3
r"""Best-effort: the W, Z, Higgs masses in the composite/walking reading (Part III).

A: m_W and m_Z computed from the condensate scale v=246 GeV and the gauge couplings --
   predicted vs observed, within ~0.4% -- plus the custodial prediction rho=1 (m_W = m_Z
   cos theta_W), the clean structural result of a doublet condensate.
B: the Higgs mass as a composite scalar -- bracketed between a light walking-dilaton
   (~0.5v) and a generic/QCD-scaled heavy sigma (~1 TeV). The observed 125 GeV sits at the
   LIGHT edge: the dilaton of the same walking the S-parameter needs.

Renders figures/pdf/electroweak_masses.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import electroweak_masses as ew

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    s = ew.summary()
    lo, hi = s["m_H"]["bracket_low"], s["m_H"]["bracket_high"]
    lines = [
        "Best-effort: W, Z, Higgs masses in the composite/walking reading",
        "=" * 64,
        f"  condensate scale v = {s['v_condensate']:.1f} GeV (GENERATED, not an inserted VEV)",
        f"  m_W: predicted {s['m_W']['predicted']:.2f} GeV  observed {s['m_W']['observed']:.2f}"
        f"  (off x{s['m_W']['off']:.3f})",
        f"  m_Z: predicted {s['m_Z']['predicted']:.2f} GeV  observed {s['m_Z']['observed']:.2f}"
        f"  (off x{s['m_Z']['off']:.3f})",
        f"  rho = m_W^2/(m_Z^2 cos^2 theta_W) = {s['rho']['predicted']:.4f}  "
        f"(custodial PREDICTION = 1; observed ~0.1%)",
        f"  m_W = m_Z cos theta_W = {s['m_W_custodial']:.2f} GeV  (the rho=1 relation)",
        "",
        f"  Higgs (composite scalar, bracketed):",
        f"    walking dilaton ~0.5v = {s['m_H']['walking_dilaton']:.0f} GeV   "
        f"(light edge)",
        f"    QCD-scaled heavy sigma ~5.4v = {hi:.0f} GeV   (excluded heavy edge)",
        f"    OBSERVED = {s['m_H']['observed']:.1f} GeV  (m_H/v = {s['m_H']['ratio_obs']:.3f})"
        f"  -> at the LIGHT/walking edge",
        "",
        "READING: m_W, m_Z come straight from the condensate scale + the gauge couplings",
        "(within 0.4%), and rho=1 is a genuine prediction of a custodial doublet condensate",
        "-- no inserted Higgs VEV, one fewer parameter. The 125 GeV scalar is COMPOSITE; its",
        "mass is a bound state (owed to lattice, L3), but its LIGHTNESS is explained: it is the",
        "pseudo-dilaton of the near-conformal 'walking' dynamics -- the same walking that",
        "pushes the electroweak S below the LEP bound. Heavy generic technicolor is excluded;",
        "the observed Higgs sits exactly where walking puts it.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "electroweak_masses.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.4))

    # ---- A: W, Z predicted vs observed ----
    labels = ["W", "Z"]
    pred = [s["m_W"]["predicted"], s["m_Z"]["predicted"]]
    obs = [s["m_W"]["observed"], s["m_Z"]["observed"]]
    x = np.arange(len(labels))
    axA.bar(x - 0.18, obs, width=0.34, color="C7", edgecolor="k", lw=0.6, label="observed")
    axA.bar(x + 0.18, pred, width=0.34, color="C0", edgecolor="k", lw=0.6,
            label=r"ED: $\frac{1}{2}g v$, $\frac{1}{2}\sqrt{g^2+g'^2}\,v$")
    for xi, p, o in zip(x, pred, obs):
        axA.text(xi + 0.18, p + 1, f"{p:.1f}", ha="center", fontsize=8.5, color="C0")
        axA.text(xi - 0.18, o + 1, f"{o:.1f}", ha="center", fontsize=8.5, color="0.3")
    axA.set_xticks(x)
    axA.set_xticklabels([r"$m_W$", r"$m_Z$"], fontsize=12)
    axA.set_ylabel("mass (GeV)")
    axA.set_ylim(0, 105)
    axA.set_title(f"W, Z from the condensate scale (within 0.4%)\n"
                  r"custodial $\rho=1$: $m_W = m_Z\cos\theta_W$", fontsize=11)
    axA.legend(fontsize=9, loc="upper left")
    axA.grid(True, axis="y", alpha=0.2)

    # ---- B: the Higgs bracket (log scale) ----
    axB.axhspan(lo, hi, color="0.85", alpha=0.7, label="composite-scalar bracket")
    axB.axhline(s["m_H"]["walking_dilaton"], color="C2", lw=1.6, ls="--",
                label=r"walking dilaton $\sim 0.5v$ (light)")
    axB.axhline(hi, color="C3", lw=1.6, ls="--",
                label=r"QCD-scaled $\sigma\sim 5.4v$ (heavy, excluded)")
    axB.axhline(s["m_H"]["observed"], color="C1", lw=2.6,
                label=f"observed Higgs = {s['m_H']['observed']:.1f} GeV")
    axB.scatter([0], [s["m_H"]["observed"]], color="C1", s=90, zorder=5, edgecolor="k")
    axB.annotate("at the light/walking edge\n(dilaton of the walking\nthe S-parameter needs)",
                 (0, s["m_H"]["observed"]), textcoords="offset points", xytext=(20, -10),
                 fontsize=8.5, color="C1")
    axB.set_yscale("log")
    axB.set_xlim(-1, 2)
    axB.set_xticks([])
    axB.set_ylabel("Higgs / composite-scalar mass (GeV)")
    axB.set_title("The Higgs as a composite scalar (bracketed)\n"
                  "observed sits at the light edge -> walking", fontsize=11)
    axB.legend(fontsize=8, loc="upper right")
    axB.grid(True, which="both", axis="y", alpha=0.15)

    fig.suptitle("W, Z, Higgs in the composite/walking reading: gauge-boson masses clean "
                 "(0.4%, $\\rho=1$); Higgs light because it walks (exact owed, L3)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(PDF_DIR, "electroweak_masses.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "electroweak_masses.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'electroweak_masses.pdf')}")


if __name__ == "__main__":
    main()
