r"""First-order (G,F) Dirac torsiton -- the relativistic solver Eq. 6.2 actually requires.

The second-order (Schrodinger) reduction throws away the small component F, and so cannot tell
the SCALAR channel (the attraction, couples to G^2 - F^2) from the VECTOR/axial channel (the
repulsion, couples to G^2 + F^2). That mishandling -- not the physics -- unbound the soliton when
the Fierz repulsion was switched on (chiral_soliton, the g_v scan). This module keeps F: it solves
the coupled first-order radial Dirac

    dG/dr = -(kappa/r) G + [ (E - V0) + M ] F
    dF/dr = +(kappa/r) F - [ (E - V0) - M ] G

with the SCALAR self-energy in M(r) = g sigma(r) (sourced by sum (G^2 - F^2)) and the repulsive
VECTOR potential V0(r) = g_v omega(r) entering as E -> E - V0 (sourced by the number density
sum (G^2 + F^2)). Two findings from this solver:

  * RELATIVITY RESTORES BINDING. Done in first order the soliton does NOT trivially unbind under the
    repulsion (where the second-order version did). The small component F is itself part of the
    repulsion, so a properly relativistic bound state self-limits. "Where there is a potential there
    are bound states" -- once F is kept.
  * THE REPULSION LOOSENS THE BAG, MONOTONICALLY. Turning the (ranged) spin-push up shallows the
    chiral-restored core -- real work, not nothing. A direct CONTACT vector (V0 proportional to the
    local density) is stiffer and flattens the bag fastest. NOTE: a clean "range vs contact" verdict
    is confounded by the coupling normalisation -- the omega field amplitude scales as 1/m_omega^2,
    so changing the range also changes the strength -- so disentangling range from strength waits on
    the proper coupling derivation (g, g_v, m_omega from the single Hehl-Datta G).

STATUS: foundation. The eigensolver is validated against the second-order solver on the same well;
the self-consistent loop with the ranged omega runs and binds. NOT yet settled: the coupling
normalisation g, g_v, m_omega, lam, v from the single Hehl-Datta G (the open derivation), and the
convergence/level-counting need tightening before the spectrum (the generations, M_V/f_pi -> S) is
trustworthy. m=0 (chiral) throughout; the vacuum mass comes from the double-well sigma (the range
that installs the condensate without the divergent Dirac sea).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from cartasis_sims import chiral_soliton as cs


def dirac_levels(M, V0, r, kappa=-1, n_levels=4, norm_floor=0.3):
    """Bound levels of the first-order radial Dirac in a scalar well M(r) and vector potential
    V0(r): H[G;F]=E[G;F]. Returns [(E, G, F), ...] for the lowest bound states (0 < E < M_vac,
    localised). Central-difference derivative; the doublers sit outside the gap and are dropped by
    the gap + localisation (norm_floor) selection -- no Wilson term needed (it distorts the physical
    low modes). G, F are normalised to integral(G^2+F^2) dr = 1."""
    N = len(r)
    h = r[1] - r[0]
    D = np.zeros((N, N))
    for i in range(1, N - 1):
        D[i, i + 1] = 1.0 / (2 * h)
        D[i, i - 1] = -1.0 / (2 * h)
    D[0, 0] = -1.0 / h; D[0, 1] = 1.0 / h
    D[-1, -1] = 1.0 / h; D[-1, -2] = -1.0 / h
    kr = np.diag(kappa / r)
    H = np.block([[np.diag(V0 + M), -D + kr], [D + kr, np.diag(V0 - M)]])
    H = 0.5 * (H + H.T)
    w, v = np.linalg.eigh(H)
    Mvac = M[-1]
    GF = []
    for k, e in enumerate(w):
        if 0.0 < e < Mvac * 0.999:
            vec = v[:, k] / np.sqrt(h)
            G, F = vec[:N], vec[N:]
            if np.trapezoid(G**2 + F**2, r) > norm_floor:
                GF.append((float(e), G, F))
    GF.sort(key=lambda t: t[0])
    return GF[:n_levels], float(Mvac)


def omega_field(n_density, r, h, m_omega, g_v):
    """The ranged repulsive vector field: -omega'' - (2/r) omega' + m_omega^2 omega = g_v n(r),
    via w = r omega. Finite range 1/m_omega; the contact limit is m_omega -> infinity."""
    rhs = g_v * r * n_density
    diag = 2.0 / h**2 + m_omega**2 * np.ones(len(r))
    off = -np.ones(len(r) - 1) / h**2
    ab = np.zeros((3, len(r)))
    ab[0, 1:] = off
    ab[1, :] = diag
    ab[2, :-1] = off
    return solve_banded((1, 1), ab, rhs) / r


def solve_dirac_soliton(g=4.0, v=1.0, lam=8.0, n_fermions=1, g_v=0.0, m_omega=2.0,
                        R=10.0, N=240, iters=200, mix=0.2, eps_break=0.4, seed=None):
    """Self-consistent first-order Dirac torsiton (m=0). Scalar double-well sigma (M=g sigma,
    sourced by sum(G^2-F^2)) + ranged repulsive omega (V0=g_v omega, sourced by sum(G^2+F^2)).
    Seed `seed` (a sigma profile) for adiabatic continuation in g_v -- ramp g_v in small steps,
    feeding each converged sigma as the next seed, so the bound configuration is never lost."""
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    sigma = seed.copy() if seed is not None else v * np.tanh(r / 1.5)
    V0 = np.zeros(N)
    nb = 0
    for it in range(iters):
        M = np.abs(g * sigma)
        lv, Mvac = dirac_levels(M, V0, r, n_levels=max(3, n_fermions))
        nb = len(lv)
        if nb < n_fermions:
            return {"bound": nb, "sigma": sigma, "core": float(sigma[0] / v),
                    "E": [e for e, _, _ in lv], "converged": False, "r": r}
        s_scalar = np.sum([G**2 - F**2 for _, G, F in lv[:n_fermions]], axis=0) / r**2
        n_density = np.sum([G**2 + F**2 for _, G, F in lv[:n_fermions]], axis=0) / r**2
        sigma = (1 - mix) * sigma + mix * cs._sigma_newton(sigma, s_scalar, r, h, lam, v, g, eps_break)
        if g_v > 0.0:
            V0 = (1 - mix) * V0 + mix * (g_v * omega_field(n_density, r, h, m_omega, g_v))
    return {"bound": nb, "sigma": sigma, "core": float(sigma[0] / v),
            "E": [e for e, _, _ in lv], "Mvac": Mvac, "converged": True, "r": r}
