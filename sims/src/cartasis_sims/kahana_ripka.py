r"""Kahana--Ripka basis for the chiral-quark-soliton: the CONVERGED torsiton sea energy and mass.

The finite-difference grid (chiral_quark_soliton) resolves the high grand-spin centrifugal barrier
badly: its sea sum crawls (~K^{-1.2}) and drifts with the box. The cure -- and the reason CQSM is
always done this way -- is to expand the quark field in free spherical-Bessel Dirac partial waves in
a box (Kahana--Ripka 1984): the barrier and the r^l origin behaviour are then EXACT, sigma.p acts
analytically through Bessel recurrences, and the vacuum subtraction is basis-consistent.

Construction (per grand-spin sector K, parity): for each j-channel (upper orbital l, lower l'=2j-l,
Dirac quantum number kappa) the radial momenta are the Dirichlet zeros of j_l, p_n = z_{l,n}/D, and
the basis states are the free massive Dirac partial waves
    g_n(r) = j_l(p_n r),   f_n(r) = -sign(kappa) (p_n/(E+M)) j_{l'}(p_n r),   E = +/- sqrt(p_n^2+M^2),
normalised on (0,D). The free Hamiltonian is diagonal (energies E); the hedgehog adds the chiral
mass deviation beta M(cos theta - 1) (same-channel radial integrals) and the hedgehog coupling
beta i gamma5 tau.rhat M sin theta (upper<->lower radial integrals weighted by the grand-spin angular
matrix tau.rhat of chiral_quark_soliton). The relative sign was fixed once, by matching the K=0
valence to the trusted hedgehog.kzero_levels.

VALIDATED (tests): free spectrum exact (+/- sqrt(p^2+M^2), parity-degenerate, no spurious in-gap);
the in-gap bound spectrum reproduces chiral_quark_soliton.ingap_levels; the sea energy CONVERGES in
grand spin (Delta ~ 1e-3 by Kmax=16, vs the grid's runaway) and is box-stable to ~2%; the small-
amplitude sea is a clean gradient term giving f_pi ~ 0.39 M (proper-time scheme, Appendix C).

RESULT. The B=1 hedgehog sea energy is positive and grows with size; the valence dives but only
saturates. Solving the model FULLY self-consistently (self_consistent_profile: the chiral angle
theta(r) determined by the quark source, not assumed) exposes a CRITICAL COUPLING. For Lambda above
~2.2 M (the weak, derived regime, f_pi/M ~ 0.4) the iteration collapses to the trivial vacuum
(theta -> 0): no stable soliton, and the B=1 sector is the unbound constituent sum ~ N_c M. For
Lambda below ~2.2 M (strong coupling) a STABLE self-consistent torsiton soliton forms -- the profile
winds theta(0) = pi -> 0 over ~2/M, the valence binds at eps_val ~ 0.7 M, and
    M_torsiton  ~=  3 N_c-th... = 3 eps_val + E_sea  ~=  3.6 M  ~=  N_c M ,
slightly above the constituent sum. So the torsiton exists as a chiral soliton when propagating
torsion (the L5 scale Lambda) is strongly enough coupled, with mass ~ N_c M; whether the physical
Lambda sits above or below the critical coupling is for the lattice (L4/L5). At this coupling only
the valence orbital binds (no excited in-gap level), so higher torsiton generations are not excited
orbitals of one soliton here -- they would be collective/radial excitations, a separate question.
"""
from __future__ import annotations

import functools

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfc, spherical_jn

from cartasis_sims import chiral_quark_soliton as cqs


@functools.lru_cache(maxsize=256)
def jl_zeros(l, n):
    """First n positive zeros of the spherical Bessel function j_l."""
    zeros, x, step = [], 0.05, 0.05
    f0 = spherical_jn(l, x)
    while len(zeros) < n:
        x1 = x + step
        f1 = spherical_jn(l, x1)
        if f0 * f1 < 0:
            zeros.append(brentq(lambda t: spherical_jn(l, t), x, x1))
        x, f0 = x1, f1
    return tuple(zeros[:n])


def _channels(K, parity):
    """(j, l_up, l_low, kappa, upper_index, lower_index) per j-channel of sector (K, parity),
    aligned to the angular-state ordering of chiral_quark_soliton.taurhat_matrix."""
    t, U, L = cqs.taurhat_matrix(K, parity)
    chans = []
    for a, (j, lu) in enumerate(U):
        llow = int(2 * j - lu)
        b = next(i for i, (jl, ll) in enumerate(L) if jl == j)   # the same-j lower (sigma.p partner)
        chans.append((j, lu, llow, cqs._kappa(j, lu), a, b))
    return chans, t


def _grid(D, Nq):
    rq = np.linspace(D / Nq, D, Nq)
    return rq, np.full(Nq, D / Nq) * rq ** 2


def _assemble(K, parity, theta, M, Nb, D, rq, w_r2):
    """Assemble sector (K, parity) and also return the basis arrays (energies E, upper/lower radial
    profiles G/F, channel index CH, the channel list, and the tau.rhat matrix t) so that densities
    can be formed from the eigenvectors. build_hamiltonian returns just the matrix."""
    chans, t = _channels(K, parity)
    mc = M * (np.cos(theta) - 1.0)
    ms = M * np.sin(theta)
    E, G, F, CH = [], [], [], []
    for ci, (j, lu, llow, kap, ai, bi) in enumerate(chans):
        s_k = -np.sign(kap)
        for pn in np.array(jl_zeros(lu, Nb)) / D:
            En = np.sqrt(pn ** 2 + M ** 2)
            gu = spherical_jn(lu, pn * rq)
            jl = spherical_jn(llow, pn * rq)
            for sign in (+1, -1):
                Es = sign * En
                f = s_k * (pn / (Es + M)) * jl
                nrm = np.sqrt(np.einsum("r,r->", gu ** 2 + f ** 2, w_r2))
                E.append(Es); G.append(gu / nrm); F.append(f / nrm); CH.append(ci)
    E = np.array(E); G = np.array(G); F = np.array(F); CH = np.array(CH)
    H = np.diag(E)
    same = CH[:, None] == CH[None, :]                              # beta M(cos-1): same channel only
    massblk = np.einsum("ar,br,r->ab", G, G, mc * w_r2) - np.einsum("ar,br,r->ab", F, F, mc * w_r2)
    H = H + np.where(same, massblk, 0.0)
    GF = np.einsum("ar,br,r->ab", G, F, ms * w_r2)                 # hedgehog: tau.rhat-weighted
    FG = np.einsum("ar,br,r->ab", F, G, ms * w_r2)
    Hh = np.zeros_like(H)
    for ci, (_, _, _, _, ai, bi) in enumerate(chans):
        for cj, (_, _, _, _, aj, bj) in enumerate(chans):
            m = (CH[:, None] == ci) & (CH[None, :] == cj)
            Hh += m * (t[bj, ai] * GF + t[bi, aj] * FG)
    H = 0.5 * ((H - Hh) + (H - Hh).T)                              # sign fixed by the K=0 match
    return H, E, G, F, CH, chans, t


def build_hamiltonian(K, parity, theta, M, Nb, D, rq, w_r2):
    """The hedgehog Dirac Hamiltonian of sector (K, parity) in the Kahana--Ripka basis (Nb momenta
    per channel, box radius D, quadrature grid rq with weights w_r2). theta is the chiral angle
    sampled on rq. Returns the dense symmetric matrix; its eigenvalues are the quark energies."""
    return _assemble(K, parity, theta, M, Nb, D, rq, w_r2)[0]


def ingap_levels(K, parity, theta_func, M=1.0, Nb=32, D=10.0, Nq=2200):
    """The discrete bound quark orbitals (|E|<M) of sector (K, parity) for the chiral angle
    theta_func(r): the valence (K=0, parity +) and the excited orbitals (candidate generations)."""
    rq, w_r2 = _grid(D, Nq)
    w = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, theta_func(rq), M, Nb, D, rq, w_r2)))
    return w[(w > -M + 1e-4) & (w < M - 1e-4)]


def _absreg(eps, Lam):
    e = np.abs(eps)
    return (Lam / np.sqrt(np.pi)) * (1.0 - np.exp(-(e / Lam) ** 2)) + e * erfc(e / Lam)


def sea_energy(theta_func, M=1.0, Lam=4.0, Nc=3, Kmax=14, Nb=32, D=10.0, Nq=2200):
    """The CONVERGED vacuum-subtracted, proper-time-regularised Dirac-sea energy of the hedgehog
    theta_func(r): -N_c sum_{K<=Kmax, parity} (2K+1) sum_{eps<0}[|eps|_reg - |eps^0|_reg]. Converges
    in Kmax (Delta ~ 1e-3 by Kmax~16) and is box-stable, unlike the finite-difference grid."""
    rq, w_r2 = _grid(D, Nq)
    th = theta_func(rq)
    vac = np.zeros_like(rq)
    tot = 0.0
    for K in range(Kmax + 1):
        for parity in (1, -1):
            ws = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, th, M, Nb, D, rq, w_r2)))
            w0 = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, vac, M, Nb, D, rq, w_r2)))
            tot += (2 * K + 1) * (_absreg(ws[ws < 0], Lam).sum() - _absreg(w0[w0 < 0], Lam).sum())
    return -Nc * tot


def soliton_mass(theta_func, M=1.0, Lam=4.0, Nc=3, Kmax=14, Nb=32, D=10.0, Nq=2200):
    """The B=1 hedgehog mass M_sol = N_c eps_val + E_sea for the chiral angle theta_func(r): the
    valence is the K=0 (parity +) in-gap level (the gap edge +M if none is bound yet). Returns
    (M_sol, eps_val, E_sea)."""
    lv = ingap_levels(0, +1, theta_func, M, Nb, D, Nq)
    ev = float(lv[np.argmin(np.abs(lv))]) if len(lv) else M
    es = sea_energy(theta_func, M, Lam, Nc, Kmax, Nb, D, Nq)
    return Nc * ev + es, ev, es


def _sector_densities(K, parity, theta, M, Lam, Nb, D, rq, w_r2):
    r"""Regularised scalar S(r) = sum psi^dag beta psi and pseudoscalar P(r) = sum psi^dag beta i
    gamma5 tau.rhat psi densities for sector (K, parity): the Dirac sea (eps<0, proper-time weight
    erfc(|eps|/Lam)) plus the K=0+ valence (weight 1). These source the chiral angle via
    theta = atan2(M_p, M_s) with (M_s, M_p) = -2 G_S (S, P)."""
    H, E, G, F, CH, chans, t = _assemble(K, parity, theta, M, Nb, D, rq, w_r2)
    w, V = np.linalg.eigh(H)
    nchan, Nq = len(chans), len(rq)

    def UW(c):
        U = np.zeros((nchan, Nq)); Wp = np.zeros((nchan, Nq))
        for ci in range(nchan):
            sel = CH == ci
            U[ci] = c[sel] @ G[sel]; Wp[ci] = c[sel] @ F[sel]
        return U, Wp

    def rho_s(U, Wp):
        return np.sum(U ** 2 - Wp ** 2, axis=0)

    def rho_p(U, Wp):
        out = np.zeros(Nq)
        for ci, (_, _, _, _, ai, bi) in enumerate(chans):
            for cj, (_, _, _, _, aj, bj) in enumerate(chans):
                out += t[bj, ai] * U[ci] * Wp[cj] + t[bi, aj] * Wp[ci] * U[cj]
        return -out                                           # match the Hamiltonian hedgehog sign
    S = np.zeros(Nq); P = np.zeros(Nq); eps_val = None
    for j in np.where(w < 0)[0]:
        U, Wp = UW(V[:, j]); wt = erfc(abs(w[j]) / Lam)
        S += wt * rho_s(U, Wp); P += wt * rho_p(U, Wp)
    if K == 0 and parity == 1:
        ingap = np.where((w > -M + 1e-4) & (w < M - 1e-4))[0]
        if len(ingap):
            jv = ingap[np.argmin(np.abs(w[ingap]))]; eps_val = float(w[jv])
            U, Wp = UW(V[:, jv]); S += rho_s(U, Wp); P += rho_p(U, Wp)
    return S, P, eps_val


def self_consistent_profile(Lam, M=1.0, Kmax=10, Nb=26, D=10.0, Nq=1600,
                            iters=28, mix=0.4, r0_init=1.2, tol=1e-3):
    r"""Solve the full chiral-quark-soliton self-consistently: the chiral angle theta(r) is NOT
    assumed but determined by the quark source, theta(r) = atan2(-P(r), -S(r)) (so theta=0 in the
    vacuum, where S<0 is the condensate and P=0), iterated to a fixed point from a B=1 seed.

    Returns (rq, theta, eps_val): the converged profile (theta(0)~pi if a soliton forms, ~0 if it
    relaxes to the vacuum) and the valence energy (None if no soliton). Across cutoffs this exposes a
    CRITICAL COUPLING: for Lam below ~2.2 M (strong coupling) a stable torsiton soliton forms; above
    it the iteration collapses to the trivial vacuum and the B=1 sector is the unbound constituent
    sum. When it binds, M_torsiton ~ N_c M ~ 3.6 M (soliton_mass on the returned profile)."""
    rq, w_r2 = _grid(D, Nq)
    theta = np.pi * np.exp(-((rq / r0_init) ** 2))
    eps_val = None
    for _ in range(iters):
        S = np.zeros(Nq); P = np.zeros(Nq); eps_val = None
        for K in range(Kmax + 1):
            for parity in (1, -1):
                s, p, e = _sector_densities(K, parity, theta, M, Lam, Nb, D, rq, w_r2)
                S += (2 * K + 1) * s; P += (2 * K + 1) * p
                if e is not None:
                    eps_val = e
        thnew = np.arctan2(-P, -S)
        thnew = np.where(thnew < -0.5, thnew + 2 * np.pi, thnew)   # stay on the B=1 winding branch
        new = (1 - mix) * theta + mix * thnew
        if np.max(np.abs(new - theta)) < tol:
            theta = new
            break
        theta = new
    return rq, theta, eps_val
