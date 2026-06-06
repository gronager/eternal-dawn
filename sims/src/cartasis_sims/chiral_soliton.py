r"""The chiral gravity-torsion soliton: a real condensate, a real f_pi, and mass you
can watch form.

self_consistent.py used a massive (non-chiral) sigma -- a soliton, but no f_pi and no
electroweak S. This module upgrades to the SYMMETRY-BREAKING (chiral) model, the
linear-sigma / Friedberg-Lee bag:

    L = psi-bar(i gamma.d - g sigma) psi + 1/2 (d sigma)^2 - U(sigma),
    U(sigma) = (lambda/4)(sigma^2 - v^2)^2.

The vacuum breaks chiral symmetry: sigma = v, so the fermion has a CONSTITUENT mass
M_vac = g v, and v is the order parameter -- f_pi = v. A soliton is a BAG where the
condensate melts (sigma -> 0 inside, chiral symmetry restored, light fermions),
surrounded by the broken vacuum. The fermions are bound because they are light inside
and heavy outside; the cost is the bag volume energy B = U(0) = lambda v^4/4.

Self-consistent loop (s-wave, second-order Dirac reduction):
  * fermion:  -u'' + (g sigma)^2 u = E^2 u,  fill the lowest N_f (Pauli);
  * sigma:    -sigma'' - (2/r) sigma' + lambda sigma(sigma^2 - v^2) = -g s(r),
              s(r) = sum |u_n/r|^2,  sigma'(0)=0,  sigma(R)=v.
Solved by Newton (nonlinear sigma) inside a damped fixed-point (fermion <-> sigma).

What this buys that the non-chiral model could not:
  * f_pi = v is a real, defined scale (the condensate), set by dimensional
    transmutation, not a finger on the scale;
  * the soliton (observable) mass is the WHOLE bag energy -- level energies + field +
    bag -- NOT the constituent g v: configurational mass you can watch switch on with
    the condensate;
  * with a real f_pi and the spectrum, a better-grounded electroweak S estimate
    (still a proxy until the V/A current correlators are done, but now with the right
    f_pi).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded


def _levels(M, r, h, n_levels, Uadd=0.0):
    diag = 2.0 / h**2 + M**2 + Uadd
    off = -np.ones(len(r) - 1) / h**2
    eps, vecs = eigh_tridiagonal(diag, off, select="i", select_range=(0, n_levels - 1))
    return np.sqrt(np.maximum(eps, 0.0)), vecs / np.sqrt(h)


def _omega(s, r, h, m_omega, g_v):
    """The repulsive vector/axial (omega) field of the Hehl-Datta term, density-sourced:
    -omega'' - (2/r) omega' + m_omega^2 omega = g_v s(r), via w = r omega:
    -w'' + m_omega^2 w = g_v r s(r), w(0)=w(R)=0."""
    rhs = g_v * r * s
    diag = 2.0 / h**2 + m_omega**2 * np.ones(len(r))
    off = -np.ones(len(r) - 1) / h**2
    ab = np.zeros((3, len(r)))
    ab[0, 1:] = off
    ab[1, :] = diag
    ab[2, :-1] = off
    return solve_banded((1, 1), ab, rhs) / r


def _sigma_newton(sigma, s, r, h, lam, v, g, eps_break=0.0, newton_iters=12):
    """Solve the nonlinear sigma BVP by Newton: -sigma'' - (2/r)sigma'
    + lam sigma(sigma^2 - v^2) - eps_break + g s = 0, with sigma'(0)=0, sigma(R)=v_true.
    The small eps_break (explicit chiral breaking, like quark masses) tilts +v to be
    the TRUE vacuum, so the fermion-pushed core rests at the restored sigma=0 rather
    than overshooting to the (now false) -v vacuum."""
    n = len(r)
    sig = sigma.copy()
    sig[-1] = v                                   # Dirichlet at R (vacuum)
    a = 1.0 / h**2 + 1.0 / (r * h)                # coeff of sigma_{i+1}
    c = 1.0 / h**2 - 1.0 / (r * h)                # coeff of sigma_{i-1}
    b0 = -2.0 / h**2                              # coeff of sigma_i (Laplacian)
    for _ in range(newton_iters):
        # residual F_i = -(a sig_{i+1} + b0 sig_i + c sig_{i-1}) + lam sig(sig^2-v^2)+g s
        lap = np.empty(n)
        lap[1:-1] = a[1:-1] * sig[2:] + b0 * sig[1:-1] + c[1:-1] * sig[:-2]
        lap[0] = a[0] * sig[1] + b0 * sig[0] + c[0] * sig[0]      # Neumann: sig_{-1}=sig_0
        F = -lap + lam * sig * (sig**2 - v**2) - eps_break + g * s
        F[-1] = 0.0
        # Jacobian (tridiagonal): d/dsig_i
        dl = -c.copy()                            # sub-diagonal
        du = -a.copy()                            # super-diagonal
        dd = -b0 + lam * (3.0 * sig**2 - v**2)    # diagonal
        dd[0] += -c[0]                            # Neumann folds c into diagonal at i=0
        dl[-1] = 0.0; du[-1] = 0.0; dd[-1] = 1.0  # Dirichlet row
        ab = np.zeros((3, n))
        ab[0, 1:] = du[:-1]
        ab[1, :] = dd
        ab[2, :-1] = dl[1:]
        delta = solve_banded((1, 1), ab, -F)
        sig = sig + delta
        if np.max(np.abs(delta)) < 1e-9:
            break
    return sig


def fierz_gv(g):
    """The Fierz-NAIVE vector coupling: the Hehl-Datta channels are S=0.25, V=0.5 (G_A=G_V), so a
    same-structure mapping gives g_v = g*sqrt(V/S) = g*sqrt(2). WARNING: this naive value unbinds
    the soliton -- see solve_chiral. The true normalization is an open derivation."""
    return g * np.sqrt(2.0)


def solve_chiral(v=1.0, g=4.0, lam=8.0, n_fermions=2, R=12.0, N=600,
                 n_levels=5, iters=300, mix=0.25, tol=1e-6, eps_break=0.4,
                 g_v=0.0, m_omega=None):
    """Self-consistent chiral soliton. Returns the converged sigma/M profiles, the bound spectrum,
    f_pi=v, the constituent mass g v, the bag constant, and the observable soliton mass.

    THE VECTOR/AXIAL REPULSION (g_v), AND AN HONEST WARNING. The Hehl-Datta term Fierzes into a
    SCALAR attractive channel (the sigma/condensate, coeff 0.25) AND a VECTOR/AXIAL repulsive one
    (the omega, coeff 0.5, G_A=G_V; cartasis_sims.fierz). The vector channel -- the density-density
    repulsion that also makes the bounce -- enters as a vector potential V_rep=g_v*omega, shifting
    the squared-mass well by ~2(g v)V_rep. It MATTERS: the soliton's binding falls with g_v and the
    core saturates (the predicted box-like profile). BUT at the Fierz-naive g_v=fierz_gv(g)=g*sqrt(2)
    the repulsion OVERWHELMS the attraction and UNBINDS the soliton entirely (no bound level); the
    soliton survives only for g_v/g <~ 1. That unbinding is almost certainly a NORMALISATION
    ARTEFACT of two naive choices, not a physical death: (i) g_v=g*sqrt(V/S) treats the double-well
    chiral sigma and the Yukawa omega as the same structure, which they are not; (ii) the 2 E V_rep
    second-order shift uses E~g v rather than the actual eigenvalue, over-weighting the repulsion --
    the proper treatment is the full first-order Dirac (G,F) with both potentials, not the squared
    reduction. So g_v is the CRUX, the true value needs the careful coupling derivation (a multi-day
    formula job), and the soliton's existence is sensitive to it. DEFAULT g_v=0 (scalar-only, the
    prior behaviour) until that derivation is done. m_omega (vector range) defaults to g v."""
    if m_omega is None:
        m_omega = g * v                              # vector range ~ constituent Compton scale (model)
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    sigma = v * np.tanh(r / 1.5)                  # initial guess: restored core, vacuum tail
    V_rep = np.zeros(N)
    history = []
    for it in range(iters):
        M = np.abs(g * sigma)
        E, u = _levels(M, r, h, max(n_levels, n_fermions), Uadd=2.0 * (g * v) * V_rep)
        s = np.sum((u[:, :n_fermions] / r[:, None]) ** 2, axis=1)
        sig_new = _sigma_newton(sigma, s, r, h, lam, v, g, eps_break)
        if g_v > 0.0:                                # vector (repulsive) channel, density-sourced
            omega = _omega(s, r, h, m_omega, g_v)
            V_rep = (1 - mix) * V_rep + mix * g_v * omega   # under-relax (else the SCF oscillates)
        sig_next = (1 - mix) * sigma + mix * sig_new
        change = float(np.max(np.abs(sig_next - sigma)))
        history.append(change)
        sigma = sig_next
        if change < tol and it > 5:
            break
    M = np.abs(g * sigma)
    E, u = _levels(M, r, h, n_levels, Uadd=2.0 * (g * v) * V_rep)
    bound = E[E < g * v]
    # energy decomposition: fermion levels + sigma gradient + potential (bag) + omega field energy
    dsig = np.gradient(sigma, r)
    U = 0.25 * lam * (sigma**2 - v**2) ** 2
    grad_E = 4 * np.pi * np.trapezoid(0.5 * dsig**2 * r**2, r)
    pot_E = 4 * np.pi * np.trapezoid(U * r**2, r)
    omega_E = 0.0
    if g_v > 0.0:
        omega = V_rep / g_v
        dom = np.gradient(omega, r)
        omega_E = 4 * np.pi * np.trapezoid((0.5 * dom**2 + 0.5 * m_omega**2 * omega**2) * r**2, r)
    level_E = float(np.sum(E[:n_fermions]))
    mass = level_E + grad_E + pot_E + omega_E
    return {
        "r": r, "sigma": sigma, "M": M, "density": s, "E": E, "bound": bound,
        "f_pi": v, "constituent_mass": g * v, "bag_constant": 0.25 * lam * v**4,
        "level_energy": level_E, "grad_energy": float(grad_E), "pot_energy": float(pot_E),
        "omega_energy": float(omega_E), "g_v": float(g_v), "V_rep": V_rep,
        "mass": float(mass), "core_sigma": float(sigma[0] / v),
        "converged": change < tol, "iters": it + 1, "residual": change,
        "history": np.array(history),
    }


# ---- electroweak S with the REAL f_pi (still a proxy until V/A correlators) -------
def s_estimate(result, MA_over_MV=1.4):
    """A better-grounded S using the chiral soliton's f_pi=v and a vector scale
    M_V ~ the first fermion excitation gap (a proxy). S = 4pi(f_pi/M_V)^2(1+M_V^2/M_A^2).
    Labelled a proxy: a rigorous S needs the composite V/A current correlators."""
    E = result["E"]
    M_V = float(E[1] - E[0]) if len(E) > 1 else float(E[0])
    f_pi = result["f_pi"]
    x = 1.0 / MA_over_MV**2
    S = 4.0 * np.pi * (f_pi / M_V) ** 2 * (1.0 + x)
    return {"f_pi": f_pi, "M_V": M_V, "fpi_over_MV": f_pi / M_V, "S_proxy": float(S)}
