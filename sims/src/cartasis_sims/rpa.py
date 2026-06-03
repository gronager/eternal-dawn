r"""RPA-dressed electroweak S: the torsion structure's advantage, honestly.

The Fierz probe (fierz.py) showed the torsion four-fermion term has G_A = G_V -- the
structural prerequisite for walking. This module resums the fermion-antifermion
bubbles (RPA) to ask what that buys for S.

WHAT IT COMPUTES. The RPA dresses the one-loop vector/axial current correlators,
Pi_X,full = Pi_X / (1 - G_X Pi_X), so the low-energy S is the leading constituent-loop
value modulated by the resonance enhancement:

    S = (N_c/6pi) * [ Pi_V'(0)/(1-G_V Pi_V(0))^2 - Pi_A'(0)/(1-G_A Pi_A(0))^2 ]
                    / ( Pi_V'(0) - Pi_A'(0) ),

anchored so S -> N_c/6pi at G -> 0 (the moments Pi(0), Pi'(0) come from the
constituent spectral functions with a cutoff Lambda; 'walking' = small chiral
breaking = large Lambda/M, where Pi_A -> Pi_V).

THE RESULT. In the walking limit the vector and axial loops equalise (Pi_A/Pi_V -> 1),
and the torsion's EQUAL couplings (G_A=G_V) keep S bounded while a QCD-like sector
(G_A < G_V) blows up: S_torsion/S_QCD falls from ~1 (strong breaking) to ~0.3 (deep
walking). So the equal couplings deliver a real, growing advantage -- the Fierz
direction is confirmed by an independent calculation.

HONEST LIMIT. This RPA is ONE-SIDED: it captures the vector-resonance enhancement
(which RAISES S -- the technicolor problem) but not the full axial-resonance /
Weinberg-sum-rule catch-up (which LOWERS S). So it OVERESTIMATES the absolute S and
cannot, by construction, show S dropping below the leading 0.16, let alone < 0.1. It
establishes the RELATIVE advantage (torsion << QCD-like, growing with walking), NOT
the absolute escape. A trustworthy absolute S < 0.1 needs the full chiral RPA (proper
sum rules, the a1 catch-up) -- likely lattice. So S stays the owed make-or-break; what
is now firmer is that the torsion structure is on the right side of it.
"""

from __future__ import annotations

import numpy as np

from . import electroweak_S as ew

_trapz = np.trapezoid


def moments(channel, M=1.0, Lam2=9.0, Nc=3, npts=60000):
    """One-loop correlator moments Pi_X(0)=int rho/s and Pi'_X(0)=int rho/s^2,
    with a UV cutoff Lambda^2 (NJL-style). channel in {'V','A'}."""
    s = np.linspace(4 * M**2 * (1 + 1e-7), Lam2, npts)
    rho = ew.spectral_V(s, M, Nc) if channel == "V" else ew.spectral_A(s, M, Nc)
    return float(_trapz(rho / s, s)), float(_trapz(rho / s**2, s))


def s_rpa(Gfac_V, Gfac_A, M=1.0, Lam2=9.0, Nc=3):
    """RPA-dressed S for vector/axial couplings G_X = Gfac_X * G_crit(V), anchored to
    the leading-order N_c/6pi. Gfac in (0,1): below the vector critical coupling."""
    PiV0, PiVp = moments("V", M, Lam2, Nc)
    PiA0, PiAp = moments("A", M, Lam2, Nc)
    Gc = 1.0 / PiV0
    GV, GA = Gfac_V * Gc, Gfac_A * Gc
    S0 = Nc / (6.0 * np.pi)
    num = PiVp / (1 - GV * PiV0) ** 2 - PiAp / (1 - GA * PiA0) ** 2
    return float(S0 * num / (PiVp - PiAp))


def compare(Lam_over_M, Gfac=0.6, qcd_axial_ratio=0.5, M=1.0, Nc=3):
    """Torsion (G_A=G_V) vs QCD-like (G_A = qcd_axial_ratio * G_V) at the same vector
    coupling. Returns S_torsion, S_qcd, their ratio, and Pi_A(0)/Pi_V(0) (walking-ness)."""
    Lam2 = (Lam_over_M * M) ** 2
    St = s_rpa(Gfac, Gfac, M, Lam2, Nc)
    Sq = s_rpa(Gfac, qcd_axial_ratio * Gfac, M, Lam2, Nc)
    PiV0, _ = moments("V", M, Lam2, Nc)
    PiA0, _ = moments("A", M, Lam2, Nc)
    return {"S_torsion": St, "S_qcd": Sq, "ratio": St / Sq,
            "PiA_over_PiV": PiA0 / PiV0, "Lam_over_M": Lam_over_M}
