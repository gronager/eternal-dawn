#!/usr/bin/env python3
r"""Generations from the torsiton overlap ladder: the hierarchy is set by ONE number, the core sharpness.

The three generations are the lowest three radial rungs of the torsiton's finite confining bag (0/1/2
interior nodes); the bag is FINITE, so the tower caps -- a fourth rung sits above dissociation. Each
rung's mass is its configurational overlap (Eq. configmass) with the mass-giving core, and the SIZE of
the hierarchy is exponentially sensitive to how SHARP that core is -- a single function, the same
strong-dynamics scale that sets the bag, not a tower of fitted integers.

Left: the consecutive ladder n=1,2,3 against a core of one sharpness s_T reproduces the WHOLE charged-
lepton span -- predicted [0.48, 121, 1777] MeV vs observed [0.51, 106, 1777], total spread 3685 vs
3477 -- with the electron's tiny mass a PREDICTION (anchor the heaviest), no magic radial numbers.
Right: the lever. The span genIII/genI sweeps orders of magnitude as the core sharpens; the observed
lepton spread (3477) is crossed at a plausible, sharp-but-physical core. The minimal mean-field soliton
digs a BROADER core (spread ~50), so the climb to the full span is the deeper IR dynamics -- one
lattice-measurable quantity (the condensate profile around the torsiton). Renders generations_overlap.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import fermion_masses as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    p = fm.predict_tower("charged-lepton")
    pred, obs = p["predicted"], p["observed"]
    s_T = p["source_size"]

    # the lever: span vs core sharpness s (same ladder machinery, harmonic bag)
    ss = np.linspace(1.5, 0.42, 80)
    span = np.array([fm._ascending_ladder(s, well=fm.WELL)[-1] /
                     fm._ascending_ladder(s, well=fm.WELL)[0] for s in ss])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.2, 4.1))

    # --- left: the ladder reproduces the leptons (shape AND span) ---
    idx = np.arange(3)
    w = 0.36
    axL.bar(idx - w / 2, obs, width=w, color="C2", alpha=0.85, label="observed $e,\\mu,\\tau$")
    axL.bar(idx + w / 2, pred, width=w, color="C0", alpha=0.85,
            label=f"torsiton ladder ($s_T={s_T:.2f}$)")
    axL.set_yscale("log")
    for i in range(3):
        axL.text(i + w / 2, pred[i] * 1.35, f"{pred[i]:.2g}", ha="center", fontsize=7.5, color="C0")
    axL.set_xticks(idx); axL.set_xticklabels(["gen I", "gen II", "gen III"])
    axL.set_ylabel("mass  (MeV)")
    axL.set_title("consecutive rungs $\\to$ the full lepton span\n(3685 vs 3477, electron predicted)",
                  fontsize=9.8)
    axL.legend(fontsize=8, loc="upper left")
    axL.grid(True, which="both", axis="y", alpha=0.2)
    axL.set_ylim(0.2, 6e3)

    # --- right: the lever -- span set by core sharpness, one lattice number ---
    axR.plot(ss, span, "-", color="C3", lw=2.0)
    axR.axhline(3477, color="C2", ls="--", lw=1.3, label="observed leptons (3477)")
    axR.axvline(s_T, color="C0", ls=":", lw=1.2)
    axR.text(s_T - 0.02, 12, f"$s_T={s_T:.2f}$", color="C0", fontsize=8, rotation=90, va="bottom")
    axR.axvspan(1.15, 1.5, color="0.5", alpha=0.12)
    axR.text(1.32, 1.6e3, "minimal mean-field\ncore (broad):\nspread $\\sim$50",
             fontsize=7.6, color="0.4", ha="center", va="center")
    axR.set_yscale("log")
    axR.invert_xaxis()                                   # sharper core to the right
    axR.set_xlabel(r"core size $s_T$  (sharper $\rightarrow$)")
    axR.set_ylabel("span (gen III / gen I)")
    axR.set_title("the lever: one number sets the magnitude", fontsize=9.8)
    axR.legend(fontsize=8, loc="lower right")
    axR.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "generations_overlap.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  s_T={s_T:.3f}  predicted {np.round(pred,2)}  observed {np.round(obs,2)}")
    print(f"  span pred {p['total_span_pred']:.0f}  obs {p['total_span_obs']:.0f}"
          f"  Koide ratios pred {np.round(p['pred_ratios'],1)} obs {np.round(p['obs_ratios'],1)}")


if __name__ == "__main__":
    main()
