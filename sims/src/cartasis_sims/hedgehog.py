r"""The hedgehog: grand-spin K=0 quark levels in a chiral soliton -- the binding that was missing.

The single-channel valence (one kappa) was marginal: it barely bound and the bag filled in
(dirac_sea.solve_sea_soliton). The chiral soliton binds the nucleon a different way -- the HEDGEHOG:
the chiral field is a position-dependent chiral rotation M_vac * exp(i gamma5 tau.rhat theta(r))
(the pion pointing radially in isospin), which locks spin and isospin into the conserved GRAND SPIN
K = J + I. In the K=0 sector a level DIVES out of the continuum into the mass gap as the chiral
angle theta(r) grows from 0 to pi -- the valence quark, deeply bound. This module computes that.

Grand-spin K=0 reduction (two radial functions, validated by the theta=0 free limit -- no in-gap
level -- and the diving):
    H = [[ M cos theta,  -d/dr - 2/r + M sin theta ],
         [ d/dr + M sin theta,  -M cos theta ]]
de-doubled with the SLAC derivative (chiral). theta=0 -> no bound level (free); theta0 -> pi -> a
level at E ~ 0.2 M (deeply bound).

STATUS: the binding MECHANISM, demonstrated. The level dives and binds where the single channel was
marginal. NOT yet done: the self-consistent hedgehog (the chiral angle theta(r) determined by the
quark sea+valence source, the soliton mass minimised) -- the full chiral-quark-soliton that would
give the torsiton mass. But the missing depth is now found.
"""
from __future__ import annotations

import numpy as np

from cartasis_sims.dirac_sea import slac_derivative


def hedgehog_profile(theta0, r, r0=1.5):
    """A chiral-angle profile theta(r) = theta0 * exp(-(r/r0)^2): theta(0)=theta0, theta(inf)=0.
    theta0 = pi is one unit of winding (baryon number 1)."""
    return theta0 * np.exp(-((np.asarray(r, dtype=float) / r0) ** 2))


def kzero_levels(theta, r, M=1.0):
    """The grand-spin K=0 quark levels INSIDE the mass gap (-M, M) for a hedgehog with chiral-angle
    profile theta(r) -- the bound (valence + diving) states. Two-channel radial Dirac, de-doubled
    (SLAC). Returns the in-gap eigenvalues, sorted."""
    theta = np.asarray(theta, dtype=float)
    N = len(r)
    h = r[1] - r[0]
    D = slac_derivative(N, h)
    c = np.diag(M * np.cos(theta))
    s = np.diag(M * np.sin(theta))
    twor = np.diag(2.0 / r)
    H = np.block([[c, -D - twor + s], [D + s, -c]])
    H = 0.5 * (H + H.T)
    w = np.linalg.eigvalsh(H)
    return np.sort(w[(w > -M) & (w < M)])


def valence_energy(theta0, r, M=1.0, r0=1.5):
    """The deepest K=0 level (the valence quark) for a hedgehog of amplitude theta0; np.nan if none
    is bound yet. Dives from +M (theta0=0) toward ~0.2 M (theta0=pi)."""
    levels = kzero_levels(hedgehog_profile(theta0, r, r0), r, M)
    if len(levels) == 0:
        return float("nan")
    return float(levels[np.argmin(np.abs(levels))])
