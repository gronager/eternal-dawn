r"""De-doubled chiral Dirac sea -- the SLAC derivative (doubling-free AND chiral).

The central-difference derivative doubles. The Wilson term that cures doubling does so by adding an
EXPLICIT MASS, which breaks chiral symmetry and distorts the m=0 spectrum -- the lattice
chiral-fermion problem (Nielsen-Ninomiya) in 1D. The SLAC derivative is doubling-free AND chiral, so
it is the right tool for the regularised Dirac-sea sum the chiral-quark-soliton needs.

STATUS: foundation. This validates the FREE-FIELD sea condensate STRUCTURE against Appendix C,
<qbar q> = -4 N_c M I0 -- negative, growing in magnitude with M. NOT yet done: the exact
regulator-matched normalisation to I0, and the vacuum-subtracted sea in the SOLITON background with
self-consistency (the chiral-quark-soliton mode sum -- the calculation that decides whether the
torsiton binds at its own couplings or hands the spectrum to the lattice).
"""
from __future__ import annotations

import numpy as np


def slac_derivative(N, h):
    """SLAC (doubling-free, chiral) first-derivative matrix on a uniform grid:
    D_ij = (-1)^(i-j)/((i-j) h) for i != j, 0 on the diagonal."""
    i = np.arange(N)
    diff = i[:, None] - i[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        D = np.where(diff != 0, ((-1.0) ** diff) / (diff * h), 0.0)
    return D


def free_sea_scalar_sum(M, R=8.0, N=160, kappas=(-1, 1, -2, 2, -3, 3), Lambda_over_M=3.0):
    """Sum of (G^2 - F^2) over the filled (negative-energy, |E| < Lambda) Dirac modes of the FREE
    field (uniform mass M), de-doubled with the SLAC derivative and summed over partial waves with
    degeneracy 2|kappa|. A STRUCTURE proxy for <qbar q> = -4 N_c M I0 (negative, growing with M); it
    is not yet normalised to I0 -- matching the box sum to the 4D-cutoff loop is the next step."""
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    D = slac_derivative(N, h)
    tot = 0.0
    for k in kappas:
        kr = np.diag(k / r)
        H = np.block([[M * np.eye(N), -D + kr], [D + kr, -M * np.eye(N)]])
        H = 0.5 * (H + H.T)
        w, v = np.linalg.eigh(H)
        mask = (w < 0) & (np.abs(w) < Lambda_over_M * M)
        for j in np.where(mask)[0]:
            vec = v[:, j]
            G, F = vec[:N], vec[N:]
            tot += 2 * abs(k) * np.sum(G**2 - F**2)
    return float(tot)
