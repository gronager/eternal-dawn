r"""The self-consistent gravity-torsion soliton (the iteration, not the toy).

soliton.py put one fermion in a FIXED well. The real object is self-consistent: the
fermions fill the levels (Pauli), their density sources the four-fermion (sigma)
field, that field is the well, and you iterate to a fixed point. This module runs
that Hartree loop and returns a genuine self-bound soliton -- the gate calculation
(Part II target #1) done properly.

Scheme (relativistic Hartree, s-wave, second-order reduction). For a scalar field
the upper radial component obeys a Schrodinger-like equation with eps = E^2 and an
effective potential M(r)^2 (kappa(kappa+1)=0 for s_{1/2}):

    -u'' + M(r)^2 u = E^2 u ,     M(r) = m0 - g sigma(r).

Fill the lowest N_f levels (Pauli). Their scalar density s(r) = sum |u_n/r|^2 sources
the sigma field through a Yukawa (gradient-regulated, range 1/m_sigma) four-fermion
interaction,

    -sigma'' - (2/r) sigma' + m_sigma^2 sigma = g s(r) ,

which localizes the well (Friedberg-Lee/Walecka-type). Iterate M -> levels -> s ->
sigma -> M with under-relaxation until the well stops moving. Output: the converged
profile, the self-bound spectrum, and the soliton mass.

Honest scope: this is a single-channel, second-order (no spin-orbit), number-density-
sourced, NON-chiral (massive-sigma, not a symmetry-breaking double well) Hartree. It
demonstrates the GATE -- a genuine self-bound soliton exists (Part II target #1),
not assumed -- but it does NOT give a reliable electroweak S: that needs the chiral
(symmetry-breaking) model and its composite vector/axial CURRENT correlators. The
f_pi/M_V fields returned here are crude diagnostics, not trustworthy S inputs.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal


def _solve_levels(M, r, h, n_levels):
    """Lowest bound levels of -u'' + M^2 u = E^2 u (Dirichlet at 0 and r_max)."""
    diag = 2.0 / h**2 + M**2
    off = -np.ones(len(r) - 1) / h**2
    eps, vecs = eigh_tridiagonal(diag, off, select="i", select_range=(0, n_levels - 1))
    # physical normalisation: integral u^2 dr = 1  ->  divide discrete vec by sqrt(h)
    u = vecs / np.sqrt(h)
    E = np.sqrt(np.maximum(eps, 0.0))
    return E, u


def _solve_sigma(s, r, h, m_sigma, g):
    """Solve -sigma'' - (2/r)sigma' + m_sigma^2 sigma = g s(r), via w = r sigma:
    -w'' + m_sigma^2 w = g r s(r),  w(0)=w(R)=0."""
    rhs = g * r * s
    diag = 2.0 / h**2 + m_sigma**2 * np.ones(len(r))
    off = -np.ones(len(r) - 1) / h**2
    from scipy.linalg import solve_banded
    ab = np.zeros((3, len(r)))
    ab[0, 1:] = off
    ab[1, :] = diag
    ab[2, :-1] = off
    w = solve_banded((1, 1), ab, rhs)
    return w / r


def solve_soliton(m0=1.0, g=4.0, m_sigma=1.0, n_fermions=1, R=9.0, N=600,
                  n_levels=5, iters=200, mix=0.3, tol=1e-6):
    """Run the self-consistent loop. Returns a dict with the converged well M(r),
    sigma(r), the filled density, the spectrum E_n, the soliton mass, and S proxies."""
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    sigma = 0.6 * np.exp(-(r / 1.2) ** 2)          # initial guess
    history = []
    for it in range(iters):
        M = m0 - g * sigma
        E, u = _solve_levels(M, r, h, max(n_levels, n_fermions))
        # fill the lowest n_fermions levels (Pauli); scalar source s(r) = sum |u_n/r|^2
        s = np.sum((u[:, :n_fermions] / r[:, None]) ** 2, axis=1)
        sigma_new = _solve_sigma(s, r, h, m_sigma, g)
        sigma_next = (1 - mix) * sigma + mix * sigma_new
        change = np.max(np.abs(sigma_next - sigma))
        history.append(change)
        sigma = sigma_next
        if change < tol:
            break
    M = m0 - g * sigma
    E, u = _solve_levels(M, r, h, n_levels)
    bound = E[E < m0]
    # soliton mass ~ sum of filled level energies + sigma field energy
    dsig = np.gradient(sigma, r)
    field_energy = 4 * np.pi * np.trapezoid(
        (0.5 * dsig**2 + 0.5 * m_sigma**2 * sigma**2) * r**2, r)
    mass = float(np.sum(E[:n_fermions]) + field_energy)
    # rough composite scales: f_pi ~ peak sigma (condensate), M_V ~ first excitation gap
    f_pi = float(np.max(sigma))
    M_V = float(E[1] - E[0]) if len(E) > 1 else float(E[0])
    return {
        "r": r, "M": M, "sigma": sigma, "density": s, "E": E, "bound": bound,
        "mass": mass, "field_energy": float(field_energy),
        "converged": change < tol, "iters": it + 1, "residual": change,
        "history": np.array(history),
        "f_pi": f_pi, "M_V": M_V, "MV_over_fpi": M_V / max(f_pi, 1e-9),
        "interior_flatness": float(M[0] / m0),   # how deep the chiral-restored core is
    }
