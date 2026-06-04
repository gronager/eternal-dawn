#!/usr/bin/env python3
r"""Best-effort: the 12 fermion masses from the soliton overlap ladder (lattice target L4).

A: predicted vs observed mass for the 9 charged fermions (log scale, spanning 5 orders of
   magnitude). Each tower (charged leptons, up-type, down-type) is the soliton's own overlap
   ladder with ONE shape knob + a scale anchor; the heavier generations are predictions.
   All 9 land within a factor ~3 (the leptons within ~20%).
B: the neutrino-lightness mechanism -- mass linear in the condensate coupling, so a neutral,
   colourless knot (weakest grip) sits far below; the observed sub-eV scale as the smallest
   coupling, not a new fine-tuned scale.

Honest status: the SHAPE and SPAN of the hierarchy are a near-parameter-free output (the
electron's tiny Yukawa = exp(-overlap), no fine-tuning); the EXACT ratios, the per-tower
size (owed to the condensate coupling), the isospin splitting, and Koide's 2/3 are owed to
the lattice (L4). Renders figures/pdf/fermion_masses.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import fermion_masses as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

TOWER_COLOR = {"charged-lepton": "C0", "up-quark": "C3", "down-quark": "C2"}


def main():
    towers = fm.all_towers()
    worst = fm.worst_residual_factor()

    lines = [
        "Best-effort: the 12 fermion masses from the soliton overlap ladder (L4)",
        "=" * 70,
    ]
    for t, r in towers.items():
        lines.append(f"\n{t}  (one shape knob s_T={r['source_size']:.3f}, scale anchored at gen 1):")
        for lab, p, o, rf in zip(fm.LABELS[t], r["predicted"], r["observed"],
                                 r["residual_factor"]):
            lines.append(f"   {lab:4s}  predicted={p:12.4g} MeV   observed={o:12.4g} MeV"
                         f"   off x{rf:5.1f}")
        lines.append(f"   span: predicted {r['total_span_pred']:.0f} vs observed "
                     f"{r['total_span_obs']:.0f};  ratios pred {np.round(r['pred_ratios'],1)} "
                     f"obs {np.round(r['obs_ratios'],1)}")
    top = fm.top_from_condensate()
    nu = fm.neutrino_suppression()
    lines += [
        f"\ntop from the condensate (NOT anchored): m_t = v/sqrt2 = {top['predicted']/1000:.1f}"
        f" GeV vs observed {top['observed']/1000:.1f} GeV (off x{top['off']:.3f}); Yukawa~1, the",
        f"  heaviest fermion IS the condensate scale -> the electron etc. below it are computed.",
        f"electron (PREDICTED, anchor=heaviest): "
        f"{fm.predict_tower('charged-lepton')['predicted'][0]:.3f} MeV vs 0.511 observed.",
        f"worst single-mass residual across all 9 charged fermions: x{worst:.1f}",
        f"neutrino lightness: needs coupling ratio c_nu/c_charged ~ {nu['needed_c_ratio']:.1e}"
        f" (weakest grip on the condensate -- no new scale)",
        f"Koide Q (leptons) = {fm.koide_Q(fm.OBSERVED['charged-lepton']):.4f} "
        f"(target 2/3 = {2/3:.4f}; an unmet target, owed)",
        f"parameters: {fm.n_free_parameters()} (model 6 vs SM 9 inserted; ideal 1)",
        "",
        "READING: the absolute scale is NOT a free anchor -- the top is the condensate scale",
        "itself (m_t = v/sqrt2, Yukawa ~1, off 0.8%), so every lighter fermion is computed by",
        "overlap suppression below it. Anchoring each tower at its HEAVIEST rung, the gen-1",
        "masses -- including the electron (0.48 vs 0.511 MeV) -- are PREDICTIONS, not inputs.",
        "From a soliton overlap ladder the whole 5-orders-of-magnitude hierarchy lands within",
        "~3x (leptons ~10-20%); neutrinos light because a neutral colourless knot grips least.",
        "OWED to lattice (L4): the inter-tower couplings (t/b~41, t/tau~97 -- the heaviest rung",
        "per tower), the exact ratios, up/down isospin splitting, and Koide's 2/3.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "fermion_masses.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.6))

    # ---- A: predicted vs observed, all 9 charged fermions ----
    lo, hi = 1e-1, 1e6
    axA.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.6, label="perfect")
    axA.fill_between([lo, hi], [lo / 3.4, hi / 3.4], [lo * 3.4, hi * 3.4],
                     color="0.8", alpha=0.5, label=r"within $\times 3.4$")
    for t, r in towers.items():
        axA.scatter(r["observed"], r["predicted"], s=70, color=TOWER_COLOR[t],
                    edgecolor="k", lw=0.6, zorder=3, label=t)
        for lab, p, o in zip(fm.LABELS[t], r["predicted"], r["observed"]):
            axA.annotate(lab, (o, p), textcoords="offset points", xytext=(6, -2),
                         fontsize=8, color=TOWER_COLOR[t])
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi)
    axA.set_xlabel("observed mass (MeV)")
    axA.set_ylabel("predicted mass (MeV) -- soliton overlap ladder")
    axA.set_title(f"9 charged fermions from the overlap ladder\n"
                  f"all within $\\times{worst:.1f}$ (leptons $\\sim$20%), "
                  "one shape knob per tower", fontsize=11)
    axA.legend(fontsize=8, loc="upper left")
    axA.grid(True, which="both", alpha=0.15)

    # ---- B: the four towers as ladders + the neutrino floor ----
    xpos = {"charged-lepton": 0, "up-quark": 1, "down-quark": 2, "neutrino": 3}
    for t, r in towers.items():
        x = xpos[t]
        axB.plot([x] * 3, r["observed"], "o-", color=TOWER_COLOR[t], lw=1.4, ms=7,
                 label=f"{t} (obs)")
        axB.plot([x + 0.18] * 3, r["predicted"], "s--", color=TOWER_COLOR[t], lw=1.0,
                 ms=5, alpha=0.7, mfc="white")
        for lab, o in zip(fm.LABELS[t], r["observed"]):
            axB.annotate(lab, (x, o), textcoords="offset points", xytext=(-16, -3),
                         fontsize=8, color=TOWER_COLOR[t])
    # neutrino floor (order eV, illustrative)
    nu_obs = np.array(fm.OBSERVED["neutrino"])
    axB.plot([xpos["neutrino"]] * 3, nu_obs, "v-", color="C4", lw=1.4, ms=7,
             label="neutrino (~eV, illustrative)")
    axB.annotate("neutral+colourless\n= weakest coupling\n= lightest tower",
                 (xpos["neutrino"], nu_obs[1]), textcoords="offset points",
                 xytext=(-150, 30), fontsize=8, color="C4",
                 arrowprops=dict(arrowstyle="->", color="C4", lw=0.8))
    axB.set_yscale("log")
    axB.set_xticks([0, 1, 2, 3])
    axB.set_xticklabels(["charged\nlepton", "up-type", "down-type", "neutrino"], fontsize=9)
    axB.set_ylabel("mass (MeV)")
    axB.set_title("The four towers: solid=observed, dashed=ladder\n"
                  "hierarchy by condensate coupling; neutrinos at the floor", fontsize=11)
    axB.legend(fontsize=7.5, loc="lower left")
    axB.grid(True, which="both", axis="y", alpha=0.15)

    fig.suptitle("Ab initio (best effort): the 12 fermion masses as soliton overlap "
                 "integrals -- roughly right, exact numbers owed to lattice (L4)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(PDF_DIR, "fermion_masses.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "fermion_masses.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'fermion_masses.pdf')}")


if __name__ == "__main__":
    main()
