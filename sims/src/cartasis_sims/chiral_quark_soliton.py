r"""Chiral-quark-soliton (CQSM) machinery for the torsiton: the grand-spin Dirac Hamiltonian in the
hedgehog background, its bound (valence + excited) spectrum, and the regularised Dirac-sea energy.

This is the apparatus that would carry Eq. 6.2 from the marginal single-channel bag (dirac_sea,
hedgehog) to the actual torsiton MASS and its excited levels (the "next generations" = excited
in-gap quark orbitals in the self-consistent chiral field). The single quark field has dynamical mass
M and couples to the chiral field U = exp(i gamma5 tau.rhat theta(r)) in the hedgehog ansatz; the
conserved quantum number is the GRAND SPIN K = J + T (J = orbital + quark spin, T = isospin), with
parity. The single-particle Hamiltonian (Dirac rep, real convention, true lower component = i x F) is

    h = [[ +M cos theta ,   sigma.p (kinetic)  +  M sin theta (tau.rhat hedgehog) ],
         [ sigma.p + M sin theta (tau.rhat)   ,            -M cos theta            ]]

reduced to coupled radial channels per (K, parity) -- 2 channels for K=0, 4 for K>=1 -- with the
SLAC (doubling-free, chiral) derivative on reduced radial functions u = r psi (flat dr measure, so
the Hamiltonian is naturally Hermitian).

WHAT IS VALIDATED (exactly, see tests):
  * the grand-spin angular matrix elements <K(j'l')|tau.rhat|K(jl)>: (tau.rhat)^2 = 1 to machine
    precision for K = 0,1,2 (tau_rhat / taurhat_matrix);
  * the FREE spectrum (theta=0): no spurious in-gap level in any (K,parity), and the two parity
    sectors are exactly degenerate, as free Dirac requires (build_hamiltonian);
  * the K=0 sector reproduces the independently-trusted hedgehog.kzero_levels valence dive to ~1e-4;
  * the small-amplitude sea energy (sea_energy, REGULAR profile vanishing at the origin) is a clean
    gradient term ~theta0^2 with NO spurious mass term (chiral symmetry respected), and gives an
    effective f_pi ~ 0.44 M, consistent with the proper-time f_pi^2 = 4 N_c M^2 I_1 (Appendix C).

WHAT IS NOT YET TRUSTWORTHY (the honest wall -- see sea_energy):
  * the ABSOLUTE sea Casimir energy of the FULL (theta0=pi) soliton, hence the torsiton mass. On this
    SLAC radial grid the grand-spin sum converges only as ~K^{-1.2} (the high-K centrifugal barrier
    is poorly represented by the global SLAC derivative), and the vacuum-subtracted E_sea drifts with
    the box size (it should be box-independent for a localised soliton). Both are uncontrolled at full
    amplitude: E_sea ~ +15-25 M and still climbing, which would give an absurd unbound mass. This is
    the known reason CQSM is done in the KAHANA-RIPKA basis (discrete spherical-Bessel box states),
    which represents the centrifugal barrier exactly, converges exponentially in K, and has a clean
    box-consistent subtraction. That basis is now built in cartasis_sims.kahana_ripka (reusing the
    angular machinery here); there the sea sum CONVERGES and is box-stable, and the B=1 mass comes out
    M_torsiton ~ N_c M ~ 3 M, confirming Appendix C from an actual mode sum. Use kahana_ripka for the
    converged energetics; this module is the readable reference for the angular couplings, the
    Hamiltonian structure, and the exact-limit validations.

The VALENCE / excited-orbital spectrum (ingap_levels) is exact for a GIVEN profile -- it is the
robust part. kahana_ripka.ingap_levels gives the same spectrum in the convergent basis.
"""
from __future__ import annotations

import numpy as np
from scipy.special import erfc

from cartasis_sims.dirac_sea import slac_derivative

# --- grand-spin angular matrix elements of tau.rhat (sympy, computed once and cached) ---------------


def _Y1_me(lp, mp, mu, l, m):
    from sympy.physics.wigner import gaunt
    return (-1) ** mp * gaunt(lp, 1, l, -mp, mu, m)


def _rhat_cart(lp, mp, l, m):
    from sympy import sqrt
    c = sqrt(4 * np.pi / 3)
    yp, ym, y0 = _Y1_me(lp, mp, 1, l, m), _Y1_me(lp, mp, -1, l, m), _Y1_me(lp, mp, 0, l, m)
    from sympy import I, sqrt as _s
    return (c * (ym - yp) / _s(2), c * I * (ym + yp) / _s(2), c * y0)


def _gs_vector(K, j, l, M):
    """The grand-spin state |K M (j l)> as a dict (ml, ms, mt) -> coefficient (sympy)."""
    from sympy import Rational as R
    from sympy.physics.wigner import clebsch_gordan
    d = {}
    for a in range(int(2 * j + 1)):
        mj = R(int(round(2 * (a - j))), 2)
        for mt in (R(1, 2), R(-1, 2)):
            cg1 = clebsch_gordan(j, R(1, 2), K, mj, mt, M)
            if cg1 == 0:
                continue
            for ml in range(-l, l + 1):
                ms = mj - ml
                if ms not in (R(1, 2), R(-1, 2)):
                    continue
                cg2 = clebsch_gordan(l, R(1, 2), j, R(ml), ms, mj)
                if cg2:
                    d[(ml, ms, mt)] = d.get((ml, ms, mt), 0) + cg1 * cg2
    return d


_TAU = None
_tau_cache: dict = {}


def tau_rhat(K, lower_state, upper_state):
    """The grand-spin matrix element <K (j' l')| tau.rhat | K (j l)>, real, M-independent (a
    grand-spin scalar). lower_state, upper_state are (j, l) tuples (sympy Rational j, int l)."""
    from sympy import Rational as R, N as Num
    global _TAU
    if _TAU is None:
        _TAU = {(R(1, 2), R(1, 2)): (0, 0, 1), (R(1, 2), R(-1, 2)): (1, -1j, 0),
                (R(-1, 2), R(1, 2)): (1, 1j, 0), (R(-1, 2), R(-1, 2)): (0, 0, -1)}
    key = (K, lower_state, upper_state)
    if key in _tau_cache:
        return _tau_cache[key]
    M = R(0) if K == 0 else K
    v1 = _gs_vector(K, lower_state[0], lower_state[1], M)
    v2 = _gs_vector(K, upper_state[0], upper_state[1], M)
    tot = 0
    for (ml1, ms1, mt1), x in v1.items():
        for (ml2, ms2, mt2), y in v2.items():
            if ms1 != ms2:
                continue
            rc = _rhat_cart(lower_state[1], ml1, upper_state[1], ml2)
            tx, ty, tz = _TAU[(mt1, mt2)]
            tot += x * y * (rc[0] * tx + rc[1] * ty + rc[2] * tz)
    val = float(complex(Num(tot, 20)).real)
    _tau_cache[key] = val
    return val


def grandspin_states(K):
    """The (j, l) angular channels for grand spin K: j = K +/- 1/2, l = j +/- 1/2 (l >= 0)."""
    from sympy import Rational as R
    js = [R(1, 2)] if K == 0 else [K - R(1, 2), K + R(1, 2)]
    out = []
    for j in js:
        for l in (int(j - R(1, 2)), int(j + R(1, 2))):
            if l >= 0:
                out.append((j, l))
    return out


def _kappa(j, l):
    from sympy import Rational as R
    return -(int(2 * j + 1) // 2) if abs(l - (j - R(1, 2))) < 1e-9 else (int(2 * j + 1) // 2)


def tau3_angular(K, jl, Kp, jlp):
    """The grand-spin angular matrix element <Kp M=0 (j' l')| tau_3 | K M=0 (j l)>: tau_3 is the
    identity on space and spin and +/-1 on the isospin projection m_t, so it is non-zero only for
    l'=l, and (being an isovector) it connects grand spin K to K' = K, K+/-1. This is the
    isorotation-cranking operator whose response gives the soliton moment of inertia. Real, and (for
    a scalar combination) M-independent; evaluated here at M=0."""
    from sympy import Rational as R, N as Num
    (j, l), (jp, lp) = jl, jlp
    if l != lp:
        return 0.0
    v = _gs_vector(K, j, l, R(0))
    vp = _gs_vector(Kp, jp, lp, R(0))
    tot = 0
    for key, a in v.items():
        if key in vp:
            _, _, mt = key
            tot += vp[key] * a * (2 * mt)
    return float(complex(Num(tot, 20)).real)


def taurhat_matrix(K, parity):
    """The tau.rhat coupling matrix t[b, a] = <lower_b| tau.rhat |upper_a> for sector (K, parity)."""
    S = grandspin_states(K)
    U = [s for s in S if (-1) ** s[1] == parity]
    L = [s for s in S if (-1) ** s[1] == -parity]
    return np.array([[tau_rhat(K, sw, su) for su in U] for sw in L]), U, L


def build_hamiltonian(K, parity, theta, r, M=1.0, D=None):
    """The reduced-radial grand-spin Dirac Hamiltonian h(K, parity) in the hedgehog background
    theta(r), on the grid r (reduced functions u = r psi, SLAC derivative, flat dr). Returns the
    dense ((nU+nL)*Nr) symmetric matrix; its eigenvalues are the quark single-particle energies."""
    Nr = len(r)
    if D is None:
        D = slac_derivative(Nr, r[1] - r[0])
    t, U, L = taurhat_matrix(K, parity)
    nu, nl, n = len(U), len(L), len(U) + len(L)
    H = np.zeros((n * Nr, n * Nr))
    Ms = np.diag(M * np.cos(theta))
    Mp = M * np.sin(theta)

    def blk(a, b):
        return (slice(a * Nr, (a + 1) * Nr), slice(b * Nr, (b + 1) * Nr))

    for a in range(nu):
        H[blk(a, a)] = Ms
    for b in range(nl):
        H[blk(nu + b, nu + b)] = -Ms
    for a, su in enumerate(U):
        ka = _kappa(*su)
        for b, sw in enumerate(L):
            if sw[0] == su[0]:                       # sigma.p connects the same-j upper/lower pair
                H[blk(nu + b, a)] += D + np.diag(ka / r)
                H[blk(a, nu + b)] += -D + np.diag(ka / r)
        for b, sw in enumerate(L):
            if abs(t[b, a]) > 1e-12:                  # hedgehog: all upper<->lower via tau.rhat
                bv = np.diag(Mp * t[b, a])
                H[blk(nu + b, a)] += bv
                H[blk(a, nu + b)] += bv
    return 0.5 * (H + H.T)


def ingap_levels(K, parity, theta, r, M=1.0, D=None):
    """The discrete bound quark orbitals (|E| < M) of sector (K, parity): the valence level (K=0,
    parity +) and the excited in-gap orbitals (the candidate 'next generations'). Exact for the
    given profile theta(r)."""
    w = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, theta, r, M, D)))
    return w[(w > -M + 1e-4) & (w < M - 1e-4)]


def _absreg(eps, Lam):
    r"""Proper-time-regularised |eps|: (Lam/sqrt(pi))(1-exp(-eps^2/Lam^2)) + |eps| erfc(|eps|/Lam).
    -> |eps| as Lam->inf, -> Lam/sqrt(pi) (bounded) for |eps|>>Lam."""
    e = np.abs(eps)
    return (Lam / np.sqrt(np.pi)) * (1.0 - np.exp(-(e / Lam) ** 2)) + e * erfc(e / Lam)


def sea_energy(theta, r, M=1.0, Lam=4.0, Nc=3, Kmax=12):
    r"""The vacuum-subtracted, proper-time-regularised Dirac-sea energy of the hedgehog theta(r):
        E_sea = -N_c sum_{K<=Kmax, parity} (2K+1) sum_{eps_n<0} [ |eps_n|_reg - |eps_n^0|_reg ]
    (levels matched by sorted index against the theta=0 vacuum).

    VALIDATED only in the SMALL-AMPLITUDE limit: for a REGULAR profile (theta(0)=0, vanishing at the
    origin so the pion field is non-singular) it is a clean gradient term ~theta0^2 reproducing
    f_pi ~ 0.44 M. At FULL amplitude (theta0=pi) it is NOT converged on this grid -- the grand-spin
    sum crawls (~K^{-1.2}) and the result drifts with the box size -- so the absolute torsiton mass
    must NOT be read off from it. Use the Kahana-Ripka basis for the converged number (module
    docstring). Provided here as the validated-in-its-limit building block, with the caveat explicit."""
    D = slac_derivative(len(r), r[1] - r[0])
    vac = np.zeros_like(theta)
    tot = 0.0
    for K in range(Kmax + 1):
        for parity in (1, -1):
            ws = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, theta, r, M, D)))
            w0 = np.sort(np.linalg.eigvalsh(build_hamiltonian(K, parity, vac, r, M, D)))
            tot += (2 * K + 1) * (_absreg(ws[ws < 0], Lam).sum() - _absreg(w0[w0 < 0], Lam).sum())
    return -Nc * tot
