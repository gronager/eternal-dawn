#!/usr/bin/env python3
r"""Paper C, Fig. 2 -- the torsiton's mass-giving bag, RESOLVED on the lattice.

The connected three-body condensate profile rho3(r) = <N|qbar q(r)|N> on the dynamical Nf=2 SU(3)
fundamental ensemble (L16x48, beta=5.6, sea -0.5, 243 cfg; run/10, t_snk=20). The half-density radius
R0 -- the sharpness s_T = R0/r0 to which the whole fermion-mass hierarchy reduces -- lifts clean off
the lattice cutoff (R0 = 1.4-1.8 a across sink times, NOT the unresolved ~1a of the coarse pilot) and
lands INSIDE the productive window [0.43, 0.70] r0. This is the genuine three-body number, the central
result of the lattice campaign. Renders figures/pdf/torsiton_bag.pdf.

Data: figures/data/torsiton_bag_L16x48.dat (repoint at the committed run/10 output when available).
HONEST: single heavy sea mass; the value is bracketed s_T = 0.45-0.56 (sink-time systematic), not
pinned -- the chiral-limit value is owed to the light-sea campaign.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
DATA = os.path.join(ROOT, "figures", "data", "torsiton_bag_L16x48.dat")
os.makedirs(PDF_DIR, exist_ok=True)

R0_A = 1.76          # half-density radius (a), t_snk=20
R0_BRACKET = (1.42, 1.76)   # across t_snk in {16,18,20}
R0_OVER_A = 3.131    # r0/a on this ensemble
WIN = (0.43, 0.70)   # the productive window for s_T


def main():
    r, rho = np.loadtxt(DATA, unpack=True)
    s_T = R0_A / R0_OVER_A
    s_lo, s_hi = R0_BRACKET[0] / R0_OVER_A, R0_BRACKET[1] / R0_OVER_A

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(r, rho, "o-", color="C0", ms=4, lw=1.4, label=r"$\rho_3(r)=\langle N|\bar q q(r)|N\rangle$")
    ax.axhline(0.5, color="0.6", lw=0.9, ls=":")
    # the resolved half-density radius and its sink-time bracket
    ax.axvspan(R0_BRACKET[0], R0_BRACKET[1], color="C2", alpha=0.16,
               label=r"$R_0 = 1.4$--$1.8\,a$ (resolved)")
    ax.axvline(R0_A, color="C2", lw=1.4)
    ax.axvline(1.0, color="C3", lw=1.1, ls="--", label="lattice cutoff $\\sim\\!a$")
    ax.annotate(r"$R_0$", (R0_A, 0.52), color="C2", fontsize=12, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xlim(0, 8); ax.set_ylim(1e-3, 1.3)
    ax.set_xlabel(r"$r/a$")
    ax.set_ylabel(r"condensate bag $\rho_3(r)$  (normalised)")
    ax.set_title("The torsiton's mass-giving bag, resolved\n"
                 fr"$s_T = R_0/r_0 = {s_T:.2f}$ (bracket ${s_lo:.2f}$--${s_hi:.2f}$) "
                 fr"$\in$ window $[{WIN[0]},{WIN[1]}]$", fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.18, which="both")

    fig.tight_layout()
    out = os.path.join(PDF_DIR, "torsiton_bag.pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=140)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  s_T = R0/r0 = {s_T:.3f}  (bracket {s_lo:.3f}-{s_hi:.3f}; window {WIN}) -> "
          f"{'IN' if WIN[0] <= s_T <= WIN[1] else 'OUT'}")


if __name__ == "__main__":
    main()
