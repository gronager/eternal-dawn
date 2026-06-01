#!/usr/bin/env python3
r"""The decisive shared-axis test: CMB axis (SCT) vs Galactic pole (systematic).

SCT predicts a single preferred axis -- our progenitor hole's Kerr axis -- shared by
the CMB large-angle anomalies and the galaxy-spin handedness. The confound is that the
reported galaxy-spin axis sits near the Galactic pole, exactly where a Milky-Way
classification systematic would put it. This is a clean two-hypothesis test: is the
galaxy-spin axis drawn around the CMB axis (SCT) or the Galactic pole (systematic)?
Left: posterior odds vs measurement precision sigma, for the current geometry. Right:
the precision needed to decide, vs how far apart the two candidate poles are.
Renders figures/pdf/shared_axis.pdf, writes sims/output/shared_axis.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import axes as ax

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

CMB = ax.lb_to_vec(260.0, 60.0)        # CMB axis of evil
GPOLE = ax.lb_to_vec(0.0, 90.0)        # Galactic pole
SEP = ax.discriminating_separation(CMB, GPOLE)     # ~30 deg


def main():
    # current Shamir-like geometry: galaxy axis ~30 deg from CMB, ~0 from G-pole
    th_sct_now, th_sys_now = 30.0, 0.0

    lines = [
        "The decisive shared-axis test: CMB axis (SCT) vs Galactic pole (systematic)",
        "=" * 70,
        f"  CMB axis of evil (l,b)=(260,60); Galactic pole (b=90); separation {SEP:.0f} deg.",
        "  Two hypotheses for the galaxy-spin axis:",
        "    H_SCT: drawn around the CMB axis (parent-spin imprint)",
        "    H_sys: drawn around the Galactic pole (Milky-Way classification bias)",
        "",
        "  CURRENT geometry (galaxy axis ~30 deg from CMB, ~0 from Galactic pole):",
        "    sigma   ln(LR)   odds(SCT:sys)   favours",
    ]
    for sig in (5, 10, 15, 20, 30):
        llr = ax.log_likelihood_ratio(th_sct_now, th_sys_now, sig)
        odds = ax.posterior_odds_sct(th_sct_now, th_sys_now, sig)
        fav = "SCT" if odds > 1 else "systematic"
        lines.append(f"    {sig:2d} deg   {llr:+5.2f}    {odds:.1e}     {fav}")
    sig_dec = ax.sigma_to_decide(SEP, 0.0, 100.0)
    lines += [
        "",
        f"  precision needed to DECIDE at 100:1 odds: sigma < {sig_dec:.1f} deg",
        "  (the poles are only ~30 deg apart, so the galaxy axis must be measured to",
        "   ~10 deg or better; current samples are far noisier and a-posteriori).",
        "",
        "READING: with today's numbers the test favours the SYSTEMATIC -- the galaxy-",
        "spin axis sits on the Galactic pole, not the CMB axis. That is NOT a",
        "refutation: the measurement is dominated by Milky-Way classification bias and",
        "selection. The framework makes a sharp, falsifiable demand: a de-confounded,",
        "all-sky galaxy-handedness map (Euclid, LSST/Rubin) measuring the spin axis to",
        "sigma < ~10 deg must land it on the CMB axis (260,60), NOT the Galactic pole.",
        "If it stays on the Galactic pole at that precision, the parent-spin prediction",
        "is dead. This is the cheapest decisive test in the program -- archival + near-",
        "future survey data, no new physics needed.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "shared_axis.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.8, 5.2))

    # Panel L: posterior odds vs sigma, for two scenarios (galaxy axis on Gpole now;
    # galaxy axis on CMB axis in a de-confounded future).
    sig = np.linspace(3, 35, 300)
    odds_now = np.array([ax.posterior_odds_sct(30.0, 0.0, s) for s in sig])
    odds_fut = np.array([ax.posterior_odds_sct(0.0, 30.0, s) for s in sig])
    axL.plot(sig, odds_now, "C3", lw=2.2,
             label="current: axis on Galactic pole")
    axL.plot(sig, odds_fut, "C0", lw=2.2,
             label="if de-confounded: axis on CMB axis")
    axL.axhline(1.0, color="0.5", lw=1.0)
    axL.axhline(100, color="0.6", ls=":", lw=1.0)
    axL.axhline(0.01, color="0.6", ls=":", lw=1.0)
    axL.text(33, 130, "decisive for SCT", fontsize=7.5, color="C0", ha="right")
    axL.text(33, 0.013, "decisive for systematic", fontsize=7.5, color="C3", ha="right")
    axL.axvline(ax.sigma_to_decide(SEP, 0.0), color="0.4", ls="--", lw=1.0)
    axL.text(ax.sigma_to_decide(SEP, 0.0) + 0.5, 1e5, "need $\\sigma<10°$\nto decide",
             fontsize=8, color="0.3")
    axL.set_yscale("log")
    axL.set_xlabel(r"galaxy-spin axis precision  $\sigma$  [deg]")
    axL.set_ylabel(r"posterior odds  (SCT : systematic)")
    axL.set_title("The test is sharp once the axis is measured well\n"
                  "(current data favour the systematic; future data decides)",
                  fontsize=11)
    axL.set_ylim(1e-9, 1e9)
    axL.legend(fontsize=8.5, loc="center right")
    axL.grid(True, which="both", alpha=0.2)

    # Panel R: precision needed to decide, vs separation between poles.
    seps = np.linspace(5, 60, 300)
    sig_need = np.array([ax.sigma_to_decide(s, 0.0, 100.0) for s in seps])
    axR.plot(seps, sig_need, "C0", lw=2.2, label="$\\sigma$ for 100:1 odds")
    axR.fill_between(seps, 0, sig_need, color="C0", alpha=0.12)
    axR.axvline(SEP, color="C3", ls="--", lw=1.2)
    axR.plot([SEP], [ax.sigma_to_decide(SEP, 0.0)], "ko", ms=6)
    axR.annotate(f"our case: poles {SEP:.0f}° apart\n$\\Rightarrow$ need $\\sigma<"
                 f"{ax.sigma_to_decide(SEP,0.0):.0f}°$", xy=(SEP, ax.sigma_to_decide(SEP, 0.0)),
                 xytext=(38, 6), fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color="0.5"))
    axR.text(33, 17, "decidable\n(measured better\nthan this)", fontsize=8, color="C0")
    axR.set_xlabel(r"separation between candidate poles  [deg]")
    axR.set_ylabel(r"galaxy-axis precision needed  $\sigma$  [deg]")
    axR.set_title("How well must we measure the axis?\n"
                  "closer poles demand sharper data", fontsize=11)
    axR.set_ylim(0, 22)
    axR.legend(fontsize=8.5, loc="upper left")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(os.path.join(PDF_DIR, "shared_axis.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "shared_axis.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'shared_axis.pdf')}")


if __name__ == "__main__":
    main()
