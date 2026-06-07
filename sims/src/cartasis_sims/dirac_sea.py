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


def _local_sea(M_arr, Lambda, r, D, h, kmax, Nc):
    """Local sea condensate density (N_c/4pi) sum_kappa 2|kappa| sum_n (G_n^2-F_n^2)/r^2 over the
    filled (E<0, |E|<Lambda) modes, de-doubled (SLAC), summed over partial waves."""
    N = len(r)
    rho = np.zeros(N)
    for k in list(range(-kmax, 0)) + list(range(1, kmax + 1)):
        kr = np.diag(k / r)
        H = np.block([[np.diag(M_arr), -D + kr], [D + kr, np.diag(-M_arr)]])
        H = 0.5 * (H + H.T)
        w, v = np.linalg.eigh(H)
        for j in np.where((w < 0) & (np.abs(w) < Lambda))[0]:
            vec = v[:, j] / np.sqrt(h)
            rho += 2 * abs(k) * (vec[:N] ** 2 - vec[N:] ** 2)
    return (Nc / (4.0 * np.pi)) * rho / r**2


def soliton_sea_condensate(M_profile, Lambda, r, kmax=30, Nc=3):
    """The VACUUM-SUBTRACTED sea condensate density <qbar q>(r) in a soliton well M_profile(r),
    minus the uniform vacuum M_vac = M_profile[-1] -- the finite, physical sea-polarisation piece
    (the divergence and the box artefacts cancel in the subtraction). The vacuum condensate is
    negative, so a POSITIVE subtracted value in the core means the sea has RESTORED chiral symmetry
    there, alongside the valence -- it reinforces the bag (the chiral-quark-soliton mechanism that
    adds binding). Converged in kmax (~30 for Lambda~4M, bag~1/M); localised to the bag."""
    M_profile = np.asarray(M_profile, dtype=float)
    N = len(r)
    h = r[1] - r[0]
    D = slac_derivative(N, h)
    sol = _local_sea(M_profile, Lambda, r, D, h, kmax, Nc)
    vac = _local_sea(M_profile[-1] * np.ones(N), Lambda, r, D, h, kmax, Nc)
    return sol - vac


def _lowest_valence(M, r, D, h):
    N = len(r)
    kr = np.diag(-1.0 / r)
    H = np.block([[np.diag(M), -D + kr], [D + kr, np.diag(-M)]])
    H = 0.5 * (H + H.T)
    w, v = np.linalg.eigh(H)
    Mvac = abs(M[-1])
    best = None
    for j, e in enumerate(w):
        if 0 < e < Mvac * 0.999:
            vec = v[:, j] / np.sqrt(h)
            G, F = vec[:N], vec[N:]
            if np.trapezoid(G**2 + F**2, r) > 0.3 and (best is None or e < best[0]):
                best = (float(e), G, F)
    return best


def solve_sea_soliton(M_vac=1.0, Lambda=4.0, R=8.0, N=140, kmax=16, Nc=3,
                      iters=60, mix=0.25, seed=None):
    """Self-consistent chiral-quark-soliton: M(r) = M_vac - 2 G_S [sea_sub(r) + valence(r)], the mass
    sourced by the FULL condensate (Dirac sea + valence), with NO by-hand double-well. G_S is fixed
    by the gap equation M_vac = -2 G_S <qbar q>_vac (continuum_condensate, same Hamiltonian scheme),
    so the only inputs are M_vac (the unit) and the cutoff Lambda. Returns the converged M(r), the
    core depth M(0)/M_vac, the valence energy, and G_S.

    FINDING (leading order, single-channel valence): at the derived QCD-like couplings the sea
    reinforces an externally imposed bag (soliton_sea_condensate gives a positive, localised
    restoration), but the bag is NOT self-sustaining -- starting from M(0)=0 the loop returns
    M_new(0) ~ 0.5 > 0, so it climbs back to the trivial vacuum M=M_vac, the only stable fixed
    point. So the leading mean-field PLUS the sea does not robustly self-bind the torsiton; the
    proper chiral-quark-soliton bookkeeping (the hedgehog ansatz that binds the QCD nucleon) or the
    non-perturbative lattice is what would. The laptop has carried Eq. 6.2 to the edge -- a marginal,
    non-self-sustaining bag -- and the spectrum honestly hands off."""
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    D = slac_derivative(N, h)
    G_S = -M_vac / (2.0 * continuum_condensate(M_vac, Lambda, Nc))
    M = seed.copy() if seed is not None else M_vac * np.tanh(r / 1.0)
    Eval = None
    for _ in range(iters):
        lb = _lowest_valence(M, r, D, h)
        val = (Nc / (4.0 * np.pi)) * (lb[1] ** 2 - lb[2] ** 2) / r**2 if lb else np.zeros(N)
        Eval = lb[0] if lb else None
        sea = soliton_sea_condensate(M, Lambda, r, kmax=kmax, Nc=Nc)
        M = (1 - mix) * M + mix * (M_vac - 2.0 * G_S * (sea + val))
    return {"r": r, "M": M, "core": float(M[0] / M_vac), "valence_E": Eval, "G_S": float(G_S)}
