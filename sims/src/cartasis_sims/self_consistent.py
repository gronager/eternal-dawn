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


def _solve_levels(M, r, h, n_levels, Uadd=0.0):
    """Lowest bound levels of -u'' + (M^2 + Uadd) u = E^2 u (Dirichlet at 0 and r_max).
    Uadd is an optional repulsive vector potential (the 2 m0 V_rep shift, see solve_soliton)."""
    diag = 2.0 / h**2 + M**2 + Uadd
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


def degeneracy_trend(n_list=(1, 2, 3, 4), g=4.0, g_v=0.0, m_sigma=1.0, R=12.0, N=800,
                     mix=0.1, iters=800):
    """Cross-check of the Pauli/degeneracy physics in the FULL radial level-filling solver
    against the Thomas--Fermi EOS shortcut (fermi_ball.py). Fill n_fermions = n_list levels
    (Pauli) and report, per filling, the mean radius, the interior flatness M(0)/m0 (>0 =
    chiral-restored flat core; <0 = deep well), the central density, and the binding per
    particle. With g_v=0 (scalar attraction only) more fermions just deepen the core; with the
    vector (torsion-repulsion) channel on (g_v>0) the core stays flat and the density saturates
    -- the box-like drop fermi_ball.py predicts. Returns arrays keyed by quantity."""
    n_list = list(n_list)
    out = {k: [] for k in ("n", "r_mean", "interior_flatness", "rho0", "binding_per_particle",
                           "converged")}
    for nf in n_list:
        o = solve_soliton(g=g, g_v=g_v, m_sigma=m_sigma, n_fermions=nf, R=R, N=N,
                          n_levels=max(8, nf + 2), mix=mix, iters=iters)
        out["n"].append(nf)
        out["r_mean"].append(o["r_mean"])
        out["interior_flatness"].append(o["interior_flatness"])
        out["rho0"].append(float(o["density"][0]))
        out["binding_per_particle"].append(o["binding_per_particle"])
        out["converged"].append(bool(o["converged"]))
    return {k: (np.array(v) if k != "converged" else v) for k, v in out.items()}


def solve_soliton(m0=1.0, g=4.0, m_sigma=1.0, n_fermions=1, R=9.0, N=600,
                  n_levels=5, iters=200, mix=0.3, tol=1e-6, g_v=0.0, m_omega=None):
    """Run the self-consistent loop. Returns a dict with the converged well M(r),
    sigma(r), the filled density, the spectrum E_n, the soliton mass, and S proxies.

    The torsion (Hehl--Datta) four-fermion term has TWO mean-field channels: the attractive
    SCALAR one (sigma, the chiral condensate -> M=m0-g*sigma, binding) and a repulsive
    VECTOR/axial one (omega, the same density-density repulsion that makes the bounce). With
    only the scalar channel (g_v=0, the default, original behaviour) more fermions just deepen
    the well. Switching on the vector channel (g_v>0) adds a density-sourced repulsion
    V_rep = g_v*omega entering the level equation as the leading shift 2*m0*V_rep; this is the
    Walecka sigma-omega mean field, which SATURATES -- a flat-interior, box-like drop. m_omega
    defaults to m_sigma."""
    if m_omega is None:
        m_omega = m_sigma
    r = np.linspace(R / N, R, N)
    h = r[1] - r[0]
    sigma = 0.6 * np.exp(-(r / 1.2) ** 2)          # initial guess
    V_rep = np.zeros(N)                            # repulsive vector potential (built up over the loop)
    history = []
    for it in range(iters):
        M = m0 - g * sigma
        E, u = _solve_levels(M, r, h, max(n_levels, n_fermions), Uadd=2.0 * m0 * V_rep)
        # fill the lowest n_fermions levels (Pauli); number density source s(r) = sum |u_n/r|^2
        s = np.sum((u[:, :n_fermions] / r[:, None]) ** 2, axis=1)
        sigma_new = _solve_sigma(s, r, h, m_sigma, g)
        if g_v > 0.0:                              # vector (repulsive) channel, same Yukawa, density-sourced
            omega = _solve_sigma(s, r, h, m_omega, g_v)
            V_rep = g_v * omega
        sigma_next = (1 - mix) * sigma + mix * sigma_new
        change = np.max(np.abs(sigma_next - sigma))
        history.append(change)
        sigma = sigma_next
        if change < tol:
            break
    M = m0 - g * sigma
    E, u = _solve_levels(M, r, h, n_levels, Uadd=2.0 * m0 * V_rep)
    bound = E[E < m0]
    # soliton mass ~ sum of filled level energies + sigma (and omega) field energy
    dsig = np.gradient(sigma, r)
    field_energy = 4 * np.pi * np.trapezoid(
        (0.5 * dsig**2 + 0.5 * m_sigma**2 * sigma**2) * r**2, r)
    if g_v > 0.0:
        omega = V_rep / g_v
        dom = np.gradient(omega, r)
        field_energy += 4 * np.pi * np.trapezoid(
            (0.5 * dom**2 + 0.5 * m_omega**2 * omega**2) * r**2, r)
    mass = float(np.sum(E[:n_fermions]) + field_energy)
    # rough composite scales: f_pi ~ peak sigma (condensate), M_V ~ first excitation gap
    f_pi = float(np.max(sigma))
    M_V = float(E[1] - E[0]) if len(E) > 1 else float(E[0])
    # filled radial probability sum_n u_n^2 = r^2 s ; mean radius and binding per particle
    wprob = r**2 * s
    r_mean = float(np.trapezoid(wprob * r, r) / max(np.trapezoid(wprob, r), 1e-300))
    binding_per_particle = float(m0 - np.sum(E[:n_fermions]) / max(n_fermions, 1))
    return {
        "r": r, "M": M, "sigma": sigma, "V_rep": V_rep, "density": s, "E": E, "bound": bound,
        "mass": mass, "field_energy": float(field_energy),
        "converged": change < tol, "iters": it + 1, "residual": change,
        "history": np.array(history),
        "f_pi": f_pi, "M_V": M_V, "MV_over_fpi": M_V / max(f_pi, 1e-9),
        "interior_flatness": float(M[0] / m0),   # how deep the chiral-restored core is
        "r_mean": r_mean, "binding_per_particle": binding_per_particle,
    }
