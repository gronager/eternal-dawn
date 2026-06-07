r"""De-doubled chiral Dirac sea -- the SLAC derivative (doubling-free AND chiral).

The central-difference derivative doubles. The Wilson term that cures doubling does so by adding an
EXPLICIT MASS, which breaks chiral symmetry and distorts the m=0 spectrum -- the lattice
chiral-fermion problem (Nielsen-Ninomiya) in 1D. The SLAC derivative is doubling-free AND chiral, so
it is the right tool for the regularised Dirac-sea sum the chiral-quark-soliton needs.

STATUS: free-field number VALIDATED. With enough partial waves (kmax >~ 30 for Lambda ~ 4M) the box
sea condensate box_sea_condensate reproduces the continuum result continuum_condensate to ~3% (the
residual is finite-box R, N) -- the de-doubled sea machinery gives the right NUMBER, in the
Hamiltonian (3-momentum / energy) cutoff scheme. (The earlier mismatch was pure partial-wave
truncation, not a missing factor; note this Hamiltonian scheme differs from the 4D-Euclidean I0 of
Appendix C -- the scheme must be kept consistent through the gap equation.) NOT yet done: the
vacuum-subtracted sea in the SOLITON background with self-consistency (the chiral-quark-soliton mode
sum -- the calculation that decides whether the torsiton binds at its own couplings).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad


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


def box_sea_condensate(M, Lambda, R=9.0, N=220, kmax=30, Nc=3):
    """The free-field chiral condensate density <qbar q> from the de-doubled box mode sum: sum over
    the filled (E<0, |E|<Lambda) modes of int(G^2-F^2)dr, with partial-wave degeneracy 2|kappa|, x
    N_c colours, divided by the box volume V = (4/3) pi R^3. The eigenvector (sum v_i^2 = 1) gives
    int(G^2-F^2)dr = sum(v_G^2 - v_F^2) directly. Converges to continuum_condensate as kmax grows
    (kmax >~ 30 for Lambda ~ 4M); the residual is finite-box (R, N)."""
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    D = slac_derivative(N, h)
    V = (4.0 / 3.0) * np.pi * R**3
    tot = 0.0
    for k in list(range(-kmax, 0)) + list(range(1, kmax + 1)):
        kr = np.diag(k / r)
        H = np.block([[M * np.eye(N), -D + kr], [D + kr, -M * np.eye(N)]])
        H = 0.5 * (H + H.T)
        w, v = np.linalg.eigh(H)
        for j in np.where((w < 0) & (np.abs(w) < Lambda))[0]:
            vec = v[:, j]
            tot += 2 * abs(k) * np.sum(vec[:N] ** 2 - vec[N:] ** 2)
    return Nc * tot / V


def continuum_condensate(M, Lambda, Nc=3):
    """The free chiral condensate in the same (3-momentum / energy) scheme as the box sum:
    <qbar q> = -(N_c M / pi^2) int_0^{p_max} p^2/sqrt(p^2+M^2) dp,  p_max = sqrt(Lambda^2 - M^2)
    (the energy cutoff |E| < Lambda). The reference box_sea_condensate must reproduce."""
    pmax = np.sqrt(max(Lambda**2 - M**2, 0.0))
    integral = quad(lambda p: p**2 / np.sqrt(p**2 + M**2), 0.0, pmax)[0]
    return -Nc * M / np.pi**2 * integral
