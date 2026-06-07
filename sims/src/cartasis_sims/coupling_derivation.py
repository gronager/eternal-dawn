r"""Symbolic (SymPy) derivation of the torsiton couplings from the single Hehl-Datta scale.

Reproduces, symbolically, what njl_couplings.py uses numerically -- and what Appendix C documents:
  * the NJL loop integrals I0, I1 (4D Euclidean, sharp cutoff Lambda);
  * the free-field (constituent) chiral condensate  <qbar q> = -4 N_c M I0  (the gap source);
  * the decay constant  f_pi^2 = 4 N_c M^2 I1;
  * the coupling chain that fixes the soliton, as scheme-independent RATIOS:
        m_sigma = 2 M          (chiral NJL)        ->  lambda = m_sigma^2/(2 v^2) = 2 g^2
        m_omega = sqrt(6) M     (the rho)           ->  g_v   = sqrt(m_omega^2/(2 v^2)) = sqrt(3) g
        G_V/G_S = (g_v^2/m_omega^2)/(g^2/m_sigma^2) = 2   (the Fierz ratio, recovered)

The naive g_v = g*sqrt(2) with m_omega = g v that unbound the soliton violated the last line.
"""
from __future__ import annotations

import sympy as sp


def loop_integrals():
    """The NJL loops I0 = int d^4k_E/(2pi)^4 1/(k^2+M^2) and I1 = (... )/(k^2+M^2)^2, sharp cutoff
    Lambda. 4D Euclidean measure: int d^4k/(2pi)^4 f(k^2) = (1/16pi^2) int_0^{Lambda^2} k^2 f dk^2.
    Returns (M, Lambda^2, I0, I1) as SymPy expressions."""
    k2, M, Lam2 = sp.symbols("k2 M Lambda2", positive=True)
    measure = lambda f: sp.simplify(sp.integrate(k2 * f, (k2, 0, Lam2)) / (16 * sp.pi**2))
    I0 = measure(1 / (k2 + M**2))
    I1 = measure(1 / (k2 + M**2) ** 2)
    return M, Lam2, I0, I1


def sea_condensate():
    """The free-field chiral condensate from the filled Dirac sea: <qbar q> = -4 N_c M I0 (Dirac
    trace 4, N_c colours). This is the gap-equation source the numerical sea sum must reproduce."""
    M, Lam2, I0, _ = loop_integrals()
    Nc = sp.symbols("N_c", positive=True)
    return sp.simplify(-4 * Nc * M * I0)


def f_pi_squared():
    """f_pi^2 = 4 N_c M^2 I1."""
    M, Lam2, _, I1 = loop_integrals()
    Nc = sp.symbols("N_c", positive=True)
    return sp.simplify(4 * Nc * M**2 * I1)


def coupling_relations():
    """The scheme-independent coupling ratios: lambda, g_v, and the Fierz consistency G_V/G_S, from
    m_sigma=2M, m_omega=sqrt(6)M, f_pi=v=M/g, M=g v."""
    g, v = sp.symbols("g v", positive=True)
    M = g * v
    m_sigma = 2 * M
    lam = sp.simplify(m_sigma**2 / (2 * v**2))          # m_sigma^2 = 2 lambda v^2
    m_omega = sp.sqrt(6) * M
    g_v = sp.sqrt(m_omega**2 / (2 * v**2))               # KSRF: m_omega^2 = 2 g_v^2 f_pi^2
    GS = g**2 / m_sigma**2
    GV = g_v**2 / m_omega**2
    return {"lambda": sp.simplify(lam), "g_v": sp.simplify(g_v),
            "GV_over_GS": sp.simplify(GV / GS), "m_sigma": m_sigma, "m_omega": m_omega}


if __name__ == "__main__":
    M, Lam2, I0, I1 = loop_integrals()
    print("I0 =", I0)
    print("I1 =", I1)
    print("<qbar q>_sea =", sea_condensate())
    print("f_pi^2 =", f_pi_squared())
    for k, val in coupling_relations().items():
        print(f"{k} = {val}")
