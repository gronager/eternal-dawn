#!/usr/bin/env python3
r"""The self-consistent torsiton soliton, and the critical coupling (Appendix C).

Left: the chiral angle theta(r) of the B=1 torsiton, solved fully self-consistently in the
Kahana--Ripka chiral-quark-soliton model (the profile is determined by the quark source, not assumed)
at a strong coupling Lambda = 2.0 M -- a clean hedgehog winding theta(0)=pi -> 0 over ~2/M. Right:
the converged central angle theta(0) versus the cutoff Lambda, showing the critical coupling
Lambda_c ~ 2.2 M: below it a stable torsiton forms (theta(0)=pi), above it the soliton dissolves to
the vacuum (theta(0)=0) and the B=1 sector is the unbound constituent sum. Renders
figures/pdf/torsiton_profile.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import kahana_ripka as kr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def main():
    # the stable self-consistent profile at strong coupling
    rq, theta, eps_val = kr.self_consistent_profile(Lam=2.0, Kmax=10, Nb=26, Nq=1600, iters=30)
    Msol, ev, es = kr.soliton_mass(lambda r: np.interp(r, rq, theta), Lam=2.0, Kmax=12, Nb=26)

    # central angle vs coupling: the critical threshold
    lams = [1.8, 2.0, 2.15, 2.3, 2.5, 3.0]
    th0 = []
    for L in lams:
        _, th, _ = kr.self_consistent_profile(Lam=L, Kmax=8, Nb=24, Nq=1300, iters=24)
        th0.append(th[0])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    axL.plot(rq, theta, color="C3", lw=2.2)
    axL.axhline(np.pi, color="0.6", ls=":", lw=1)
    axL.set_xlim(0, 6)
    axL.set_ylim(-0.1, np.pi + 0.2)
    axL.set_yticks([0, np.pi / 2, np.pi])
    axL.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
    axL.set_xlabel(r"$r$  $[1/M]$")
    axL.set_ylabel(r"chiral angle  $\theta(r)$")
    axL.set_title(rf"self-consistent torsiton ($\Lambda=2M$): $M\!\simeq\!{Msol:.1f}\,M$", fontsize=10.5)
    axL.text(2.6, 2.3, rf"$\varepsilon_{{\rm val}}={ev:.2f}\,M$" + "\n" + rf"$E_{{\rm sea}}={es:.2f}\,M$",
             fontsize=9, color="0.3")
    axL.grid(True, alpha=0.2)

    axR.plot(lams, th0, "o-", color="C0", lw=1.8, ms=6)
    axR.axvspan(1.5, 2.2, color="C3", alpha=0.08)
    axR.axvspan(2.2, 3.2, color="0.5", alpha=0.08)
    axR.axvline(2.2, color="0.5", ls="--", lw=1)
    axR.set_xlim(1.7, 3.1)
    axR.set_ylim(-0.2, np.pi + 0.3)
    axR.set_yticks([0, np.pi / 2, np.pi])
    axR.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
    axR.set_xlabel(r"cutoff  $\Lambda$  $[M]$  (L5 scale)")
    axR.set_ylabel(r"central angle  $\theta(0)$")
    axR.set_title(r"critical coupling $\Lambda_c\!\simeq\!2.2\,M$", fontsize=10.5)
    axR.text(1.82, 0.35, "stable\ntorsiton", fontsize=8.5, color="C3")
    axR.text(2.5, 1.4, "no soliton\n(constituent sum)", fontsize=8.5, color="0.4")
    axR.grid(True, alpha=0.2)

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "torsiton_profile.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"wrote {out}  (M_sol={Msol:.3f}, eps_val={ev:.3f}, E_sea={es:.3f})")
    print(f"theta(0) vs Lambda={lams}: {[round(x,2) for x in th0]}")


if __name__ == "__main__":
    main()
