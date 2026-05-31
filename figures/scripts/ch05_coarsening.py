#!/usr/bin/env python3
r"""Radiative coarsening of the packed foam: Hawking siphoning between OGUs.

Once accretion halts, each OGU sits at marginal Hawking balance (bath ~ neighbours'
radiation ~ T_dS ~ 1e-30 K). Negative heat capacity makes it unstable: big (cold)
OGUs absorb net radiation from small (hot) ones, which evaporate into them. Left:
OGU masses evolving -- big grow, small die -- on the Hawking timescale. Right: the
instability, dM/dtau vs M/M_bar (evaporate below, grow above, unstable balance at 1).
Renders figures/pdf/coarsening.pdf, writes sims/output/coarsening.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import void_foam as vf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> None:
    rng = np.random.default_rng(1)
    m0 = rng.lognormal(0.0, 0.22, 14)
    m0 /= m0.mean()
    tau, M = vf.simulate_coarsening(m0, tau_max=4.5)
    alive = (M[:, -1] > 0.05).sum()

    lines = [
        "Radiative coarsening of the packed foam (Hawking siphoning)",
        "=" * 58,
        "  packed foam: no fresh void -> each OGU bathed in neighbours' Hawking",
        "  radiation (~ T_dS ~ 1e-30 K) -> M ~ M_crit, marginal balance.",
        "  NEGATIVE heat capacity (T_H ~ 1/M) makes it unstable:",
        "    M > M_bar -> colder -> absorbs net -> grows",
        "    M < M_bar -> hotter -> emits net -> evaporates (runaway)",
        "",
        f"  Monte-Carlo: {len(m0)} OGUs, initial spread std={m0.std():.2f}",
        f"  after {tau[-1]:.1f} Hawking times: {alive} survive, biggest "
        f"x{M[:, -1].max():.2f}, total mass conserved to "
        f"{M[:, -1].sum()/M[:, 0].sum():.3f}",
        "",
        "READING: yes -- once growth halts, Hawking DOES take over, exactly as the",
        "intuition says. But not as uniform death: it is redistribution. The bigger",
        "(colder) OGUs siphon the smaller (hotter) ones, which evaporate into them.",
        "The foam coarsens -- fewer, larger universes -- conserving total mass, on",
        "the Hawking timescale (~1e143 s for our-mass OGUs). So OGUs are immortal on",
        "any ordinary clock but not eternal: a sub-average OGU is, very slowly, eaten",
        "by its big neighbour. (The same coarsening as mergers, mediated by radiation",
        "instead of contact.)",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "coarsening.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # Panel L: mass trajectories.
    for i in range(len(m0)):
        col = "C0" if M[i, -1] > 0.05 else "C3"
        axL.plot(tau, M[i], col, lw=1.4, alpha=0.8)
    axL.axhline(1.0, color="0.6", ls=":", lw=0.8)
    axL.text(0.1, 1.04, "initial mean", fontsize=8, color="0.4")
    axL.text(3.0, M[:, -1].max() * 0.9, "grow\n(siphon neighbours)", fontsize=8,
             color="C0")
    axL.text(2.4, 0.06, "evaporate (die)", fontsize=8, color="C3")
    axL.set_xlabel(r"time  $\tau$  [Hawking times of the mean] $\;(\sim10^{143}\,$s$)$")
    axL.set_ylabel(r"OGU mass  $M_i / \langle M\rangle_0$")
    axL.set_title("The packed foam coarsens by Hawking siphoning\n"
                  "big (cold) grow, small (hot) evaporate into them", fontsize=11)
    axL.grid(True, alpha=0.2)

    # Panel R: the instability.
    x = np.linspace(0.5, 1.6, 300)
    dMdt = x**-2 * (x**4 - 1.0)
    axR.plot(x, dMdt, "C0", lw=2.2)
    axR.axhline(0.0, color="0.6", ls=":", lw=0.8)
    axR.axvline(1.0, color="C3", ls="--", lw=1.2)
    axR.fill_between(x, 0, dMdt, where=(x > 1), color="C0", alpha=0.15)
    axR.fill_between(x, dMdt, 0, where=(x < 1), color="C3", alpha=0.15)
    axR.text(1.28, 1.1, "grow\n(absorb)", fontsize=9, color="C0")
    axR.text(0.6, -1.6, "evaporate\n(emit)", fontsize=9, color="C3")
    axR.annotate("unstable balance\n$M=M_{\\rm bar}$ (neg. heat capacity)",
                 xy=(1.0, 0.0), xytext=(0.55, 1.6), fontsize=8, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3"))
    axR.set_xlabel(r"OGU mass  $M / M_{\rm bar}$")
    axR.set_ylabel(r"net growth  $\mathrm{d}M/\mathrm{d}\tau \propto (M/M_{\rm bar})^4-1$")
    axR.set_title("Why it siphons: negative heat capacity\n"
                  "the balance at $M_{\\rm bar}$ is unstable", fontsize=11)
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "coarsening.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
