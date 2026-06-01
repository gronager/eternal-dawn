#!/usr/bin/env python3
r"""What is true and observable: eta is doubly anchored; spin vs chirality.

Left: eta is fixed two INDEPENDENT ways that agree to ~1% -- the CMB acoustic peak
ratio (omega_b) and BBN deuterium (nuclear physics at z~1e9). The "CMB = parent
Hawking" story does not touch BBN, so eta ~ 6.1e-10 is rock solid. Right: the
spin/chirality truth table -- chirality SIGN is inherited (pure lineages), but SPIN is
local (the progenitor hole's Kerr spin), so spin_us != spin_OGU; the observable from
our spin is a preferred AXIS. Renders figures/pdf/spin_eta.pdf, writes
sims/output/spin_eta.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import spin_eta as se

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    e_cmb, e_bbn, diff = se.cmb_bbn_concordance()

    lines = [
        "What is true and observable: eta anchoring, spin vs chirality",
        "=" * 62,
        "ETA IS DOUBLY ANCHORED (not up in the air):",
        f"  CMB  (omega_b -> eta): eta10 = {e_cmb:.2f}",
        f"  BBN  (deuterium D/H) : eta10 = {e_bbn:.2f}",
        f"  agree to {diff*100:.1f}% -- two independent probes (recombination z~1100;",
        "  nucleosynthesis z~1e9). 'CMB = parent Hawking' concerns only the plasma's",
        "  ORIGIN at the bounce surface; it does NOT replace recombination and leaves",
        "  BBN untouched. In SCT eta is the parent debris dragged through the extrusion",
        "  (Ch.4), then read off identically by both probes.",
        "",
        "SPIN vs CHIRALITY (do not conflate):",
        "  INHERITED down a lineage:  chirality SIGN (matter/antimatter) -> pure lineages",
        "  LOCAL, re-drawn each gen:  spin magnitude AND axis (the progenitor hole's Kerr)",
        "",
        "  So, grading the spin claims:",
        "   * 'OGU spin lies in a narrow band'    -> TRUE as a derived SHAPE (viable",
        "     seeds cluster just above survival threshold); not pinned in absolute units.",
        "   * 'a BHU's spin follows its parent'   -> FALSE (that is chirality, not spin);",
        "     a BHU's spin is its own progenitor hole's Kerr spin, set locally.",
        "   * 'spin_us == spin_OGU'               -> FALSE; our spin is our progenitor",
        "     hole's (possibly near-extremal), not the OGU's low-band seed spin.",
        "  OBSERVABLE from our spin: a PREFERRED AXIS (progenitor hole's Kerr axis) --",
        "  the CMB 'axis of evil' and galaxy-spin handedness should share it (Ch.8).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "spin_eta.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    # Panel L: D/H vs eta10, with CMB and BBN determinations crossing.
    e = np.linspace(4.0, 9.0, 300)
    axL.plot(e, se.deuterium_abundance(e), "C0", lw=2.2,
             label=r"BBN: $10^5\,$D/H $=2.6\,(6/\eta_{10})^{1.6}$")
    obs = 2.53
    axL.axhspan(obs - 0.03, obs + 0.03, color="C2", alpha=0.25)
    axL.axhline(obs, color="C2", lw=1.3, label="observed D/H (Cooke 2018)")
    axL.axvline(e_cmb, color="C3", lw=1.6, ls="--",
                label=fr"CMB $\omega_b\Rightarrow\eta_{{10}}={e_cmb:.2f}$")
    axL.plot([e_bbn], [obs], "ko", ms=7, zorder=5)
    axL.annotate(f"BBN D/H $\\Rightarrow\\eta_{{10}}={e_bbn:.2f}$\n"
                 f"(agree to {diff*100:.0f}%)", xy=(e_bbn, obs),
                 xytext=(6.6, 3.4), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axL.set_xlabel(r"baryon-to-photon ratio  $\eta_{10}=10^{10}\,\eta$")
    axL.set_ylabel(r"primordial deuterium  $10^5\,$D/H")
    axL.set_title(r"$\eta$ is anchored two ways that agree to $\sim$1%"
                  "\nCMB (recombination) and BBN (nucleosynthesis)", fontsize=11)
    axL.set_xlim(4, 9)
    axL.legend(fontsize=8.5, loc="upper right")
    axL.grid(True, alpha=0.2)

    # Panel R: spin vs chirality truth table.
    axR.axis("off")
    axR.set_title("Spin is local; chirality is inherited\n"
                  "(do not conflate them)", fontsize=11)
    rows = [
        ("quantity", "inherited?", "set by", True),
        ("chirality sign\n(matter/antimatter)", "YES", "the OG seed's\nvorticity sign", False),
        ("spin magnitude", "no", "progenitor hole's\nlocal collapse", False),
        ("spin axis", "no", "progenitor hole's\nKerr axis", False),
        ("$\\eta$ magnitude", "yes (debris)", "parent meal,\nextrusion (Ch.4)", False),
    ]
    y = 0.92
    for label, inh, by, hdr in rows:
        wt = "bold" if hdr else "normal"
        col = "0.0" if hdr else ("C2" if inh.lower().startswith(("yes", "y")) else "C3")
        axR.text(0.02, y, label, fontsize=9, fontweight=wt, va="top")
        axR.text(0.42, y, inh, fontsize=9, fontweight=wt, va="top",
                 color="0.0" if hdr else col)
        axR.text(0.66, y, by, fontsize=8.5, fontweight=wt, va="top")
        if hdr:
            axR.plot([0.0, 1.0], [y - 0.05, y - 0.05], color="0.5", lw=0.8)
        y -= 0.175
    axR.text(0.02, 0.04,
             "Observable from our spin: a preferred AXIS\n"
             "(CMB axis-of-evil = galaxy-spin handedness)",
             fontsize=9, style="italic", color="C0", va="bottom")
    axR.set_xlim(0, 1)
    axR.set_ylim(0, 1)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "spin_eta.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "spin_eta.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'spin_eta.pdf')}")


if __name__ == "__main__":
    main()
