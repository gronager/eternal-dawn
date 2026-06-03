#!/usr/bin/env python3
r"""Parameter count: the Standard Model vs the Eternal Dawn particle sector.

The claim, conditional on the programme delivering ("banking on S" + the lattice targets):
because the forces and masses are tied to the gravity-torsion field, the SM's free
parameters become DERIVED rather than inserted. This figure makes the count visual and
honest.

A: the SM ledger -- ~19 parameters (no neutrino mass) to ~26 (with neutrinos), every block
   inserted by hand (the Yukawa sprinkling, the gauge couplings, theta_QCD).
B: the ED-particle ledger -- the SAME blocks, each traced to G,hbar,c plus ONE
   dynamically-generated scale v=f_pi (dimensional transmutation, not a free input);
   ZERO new fundamental parameters in the ideal case, with the honest residuals flagged.

Renders figures/pdf/parameter_count.pdf.
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from cartasis_sims import parameter_count as pc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# status -> colour for the ED ledger
STATUS_COLOR = {
    "derived": "C0",
    "forced": "C2",
    "owed": "C3",
}


def _status_of(text):
    if text.startswith("owed"):
        return "owed"
    if "forced" in text.split(":")[0]:
        return "forced"
    return "derived"


def main():
    sm_full = pc.sm_count(with_neutrinos=True)
    sm_min = pc.sm_count(with_neutrinos=False)
    replaced = pc.replaced_count()
    new_ideal = pc.ed_new_free_parameters_ideal()

    lines = [
        "Parameter count: the Standard Model vs the Eternal Dawn particle sector",
        "=" * 70,
        f"  Standard Model free parameters: {sm_min} (no nu mass) .. {sm_full} (with nu)",
        "    " + ", ".join(f"{n} [{c}]" for n, c in pc.SM_PARAMETERS),
        "",
        f"  ED-particle base constants (shared with all of physics, NOT new): "
        f"{', '.join(pc.ed_base_constants())}",
        f"  + ONE dynamically-generated scale v=f_pi (dimensional transmutation, not free)",
        f"  ED-particle NEW fundamental free parameters (ideal case): {new_ideal}",
        f"  SM parameters the programme aims to DERIVE: {replaced}",
        "",
        "  How each SM block is accounted for:",
    ]
    for name, n, status in pc.ED_ACCOUNTING:
        lines.append(f"    [{n}] {name}: {status}")
    lines += [
        "",
        "  Honest residuals (the gap between promise and proof):",
    ]
    for r in pc.ed_honest_residuals():
        lines.append(f"    - {r}")
    lines += [
        "",
        "READING: the structure replaces ~19-26 inserted parameters with ZERO new",
        "fundamental ones -- every mass, mixing and coupling becomes a consequence of",
        "G,hbar,c and a condensate scale the theory generates for itself. That is a",
        "dramatic reduction IN STRUCTURE. The load-bearing caveat: the NUMBERS (S<0.1,",
        "exact masses, confinement) are owed to lattice -- the count is a promise the",
        "structure supports, not yet a proof.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "parameter_count.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.6))

    # ---- A: the SM ledger as a stacked bar ----
    names = [n for n, _ in pc.SM_PARAMETERS]
    counts = [c for _, c in pc.SM_PARAMETERS]
    bottom = 0
    palette = ["C7", "C9", "C1", "C4", "C6", "C5"]
    for name, c, col in zip(names, counts, palette):
        axA.bar(0, c, bottom=bottom, width=0.55, color=col, edgecolor="k", lw=0.6)
        axA.text(0.32, bottom + c / 2, f"{name}  [{c}]", va="center", ha="left", fontsize=8.5)
        bottom += c
    axA.axhline(sm_min, color="k", ls="--", lw=1.0)
    axA.text(-0.34, sm_min + 0.2, f"{sm_min} (no $\\nu$)", fontsize=8.5, ha="left")
    axA.text(-0.34, sm_full + 0.2, f"{sm_full} (with $\\nu$)", fontsize=8.5, ha="left")
    axA.set_xlim(-0.4, 1.7)
    axA.set_ylim(0, sm_full + 2)
    axA.set_xticks([])
    axA.set_ylabel("number of free parameters inserted by hand")
    axA.set_title(f"Standard Model: {sm_min}–{sm_full} free parameters\n"
                  "(every block put in by hand)", fontsize=11)

    # ---- B: the ED accounting, block by block ----
    y = 0
    yticks, ylabels = [], []
    for name, n, status_text in reversed(pc.ED_ACCOUNTING):
        st = _status_of(status_text)
        axB.barh(y, n, height=0.62, color=STATUS_COLOR[st], edgecolor="k", lw=0.6, alpha=0.85)
        axB.text(0.15, y, name, va="center", ha="left", fontsize=9, color="white",
                 fontweight="bold")
        yticks.append(y)
        ylabels.append(f"[{n}]")
        y += 1
    axB.set_yticks(yticks)
    axB.set_yticklabels(ylabels, fontsize=8.5)
    axB.set_xlabel("SM parameters in this block (now DERIVED, not inserted)")
    axB.set_xlim(0, max(c for _, c, _ in pc.ED_ACCOUNTING) + 1)
    axB.set_title(f"Eternal Dawn particle sector: {new_ideal} new fundamental parameters\n"
                  f"(all {replaced} traced to $G,\\hbar,c$ + generated $v=f_\\pi$)",
                  fontsize=11)
    legend = [
        Patch(facecolor=STATUS_COLOR["derived"], edgecolor="k",
              label="derived (structure shown, number owed)"),
        Patch(facecolor=STATUS_COLOR["forced"], edgecolor="k",
              label="forced (fixed by the label/algebra)"),
        Patch(facecolor=STATUS_COLOR["owed"], edgecolor="k",
              label="owed (honest residual)"),
    ]
    axB.legend(handles=legend, fontsize=8, loc="lower right")
    axB.grid(True, axis="x", alpha=0.2)

    fig.suptitle("Far fewer parameters than the SM: ~19–26 inserted "
                 "→ 0 new fundamental (structure supports it; numbers owed to lattice)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(PDF_DIR, "parameter_count.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "parameter_count.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'parameter_count.pdf')}")


if __name__ == "__main__":
    main()
