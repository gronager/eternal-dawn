r"""The torsiton couplings, DERIVED from the single Hehl-Datta scale by NJL bosonization.

Eq. 6.2's four-fermion term, Fierzed into a scalar channel (G_S = 0.25 G_HD, attractive) and a
vector/axial channel (G_V = 0.5 G_HD, repulsive; cartasis_sims.fierz), is bosonized: the auxiliary
sigma and omega fields acquire their potential, masses, and kinetic terms from the fermion loop.
At leading order (m=0, chiral) the loop gives the standard NJL relations -- and they fix every
soliton coupling as a ratio, with ONE dimensionful input left, the cutoff Lambda (= the framework's
scale, lattice target L5):

  * gap equation:      fixes the constituent mass M against Lambda;
  * decay constant:    f_pi^2 = 4 N_c M^2 I_1(Lambda);   g == M / f_pi  (Goldberger-Treiman);
  * scalar:            m_sigma = 2 M  (chiral NJL),  lambda = m_sigma^2/(2 v^2) = 2 g^2;
  * vector:            m_omega = sqrt(6) M  (the rho ~ 770 MeV),  g_v = sqrt(3) g  (KSRF);
  * consistency:       G_V/G_S = (g_v^2/m_omega^2)/(g^2/m_sigma^2) = 2  (the Fierz ratio) -- exact.

So the naive g_v = g*sqrt(2) with m_omega = g v that unbound the soliton was wrong on BOTH: the
right coupling is stronger (sqrt(3) g) but the right range is much shorter (m_omega = sqrt(6) M),
and they recombine to the Fierz 2 in integrated strength -- leaving the SHAPE (and the relativistic
G^2 +/- F^2 balance) to decide binding. These are leading-order, regularisation-scheme-dependent at
the coefficient level (the exact m_sigma/M, m_omega/M can shift beyond leading order / with the
cutoff scheme); the structural results -- m_omega > m_sigma, the Fierz ratio 2, lambda ~ g^2 -- are
robust. Feed njl_couplings(...) straight into dirac_soliton.solve_dirac_soliton.
"""
from __future__ import annotations

import numpy as np


def _I1(Lambda_over_M):
    """NJL loop I_1 = int d^4k_E/(2pi)^4 1/(k^2+M^2)^2 = (1/16pi^2)[ln(1+x) - x/(1+x)], x=(L/M)^2."""
    x = float(Lambda_over_M) ** 2
    return (1.0 / (16.0 * np.pi**2)) * (np.log(1.0 + x) - x / (1.0 + x))


def njl_couplings(g=None, Lambda_over_M=None, Nc=3):
    """The consistent torsiton coupling set, in units of the constituent mass M=1. Give EITHER
    g = M/f_pi directly, OR Lambda_over_M (the cutoff in units of M) to derive g from the loop.
    Returns M, f_pi=v, g, g_v, m_sigma, m_omega, lambda, the channel strengths G_S,G_V and their
    ratio (must be 2). The overall scale is set later by fixing v to the physical value."""
    if g is None:
        if Lambda_over_M is None:
            raise ValueError("give either g or Lambda_over_M")
        fpi_over_M = np.sqrt(4.0 * Nc * _I1(Lambda_over_M))
        g = 1.0 / fpi_over_M
    g = float(g)
    M = 1.0
    v = M / g                                   # f_pi = v
    m_sigma = 2.0 * M                           # chiral NJL
    m_omega = np.sqrt(6.0) * M                  # NJL vector (rho ~ sqrt6 M)
    g_v = np.sqrt(3.0) * g                      # KSRF: g_rho = sqrt3 g
    lam = 2.0 * g**2                            # m_sigma^2 = 2 lam v^2  ->  lam = 2 g^2
    G_S = g**2 / m_sigma**2
    G_V = g_v**2 / m_omega**2
    return {"M": M, "f_pi": v, "v": v, "g": g, "g_v": g_v,
            "m_sigma": m_sigma, "m_omega": m_omega, "lam": lam,
            "G_S": G_S, "G_V": G_V, "GV_over_GS": G_V / G_S,
            "Lambda_over_M": (None if Lambda_over_M is None else float(Lambda_over_M))}
