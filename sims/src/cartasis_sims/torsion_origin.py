r"""Upstream of the NJL chain: where does the four-fermion coupling come from, is it strong enough?

`coupling_derivation.py` takes the dynamical mass M as given and derives the soliton couplings as
ratios. This module asks the prior question the *binding* bet hinges on (dirac_woods_saxon: three
generations need m_vac*r0 ~ 10): what FIXES the four-fermion coupling G_S, and does the NJL gap even
turn on? Two candidate origins for G_S differ by ~40 orders of magnitude:

  (A) GEOMETRIC / gravitational (standard Einstein-Cartan). Eliminating the non-propagating ECSK
      torsion leaves the Hehl-Datta contact term with coupling  G_S = 3 kappa/16,  kappa = 8 pi G
      = 1/M_Pl^2.  This is Newton's G -- Planck-suppressed.

  (B) INDEPENDENT / strong (the framework's bet -- 'torsion is its own thing'). Torsion is a strong
      IR field with its own coupling  G_S = c/Lambda^2,  Lambda the substrate scale (the same one
      that sets sigma = 2 pi v^2), c of order one or larger.

The scalar NJL gap M = 8 G_S N_c M I0 has a non-trivial root iff  G_S > G_crit = 2 pi^2/(N_c Lambda^2).
We compute G_S/G_crit for both origins and the coupling strength m_vac*r0 each implies.

VERDICT (printed by report()): the soliton -- hence ALL infrared matter in this framework -- requires
(B). The geometric coupling (A) is sub-critical by ~(Lambda/M_Pl)^2 ~ 1e-37, so it gives M = 0 exactly
(no condensate, no torsiton). So the framework's matter is NOT a consequence of Einstein-Cartan torsion;
it REQUIRES torsion to be a strong field decoupled from G -- a central dynamical postulate. This is the
algebraic version of the caution 'spin/torsion is its own thing.'
"""
from __future__ import annotations

import sympy as sp

from . import coupling_derivation as cd


def njl_critical_coupling():
    """G_crit for the scalar NJL gap: 1 = 8 G N_c I0(M->0), I0(0)=Lambda^2/(16 pi^2)
    -> G_crit = 2 pi^2/(N_c Lambda^2). Returns (G_crit(symbolic in N_c, Lambda2), I0_at_0)."""
    M, Lam2, I0, _ = cd.loop_integrals()
    I0_0 = sp.simplify(sp.limit(I0, M, 0))            # = Lambda2/(16 pi^2)
    Nc = sp.symbols("N_c", positive=True)
    G_crit = sp.simplify(1 / (8 * Nc * I0_0))
    return G_crit, I0_0


def fork():
    """G_S/G_crit for the two origins, symbolically."""
    Nc, Lam, MPl, c = sp.symbols("N_c Lambda M_Pl c", positive=True)
    Lam2 = Lam**2
    G_crit = 2 * sp.pi**2 / (Nc * Lam2)               # from njl_critical_coupling (Lambda2 -> Lam^2)
    G_grav = sp.Rational(3, 16) / MPl**2              # Hehl-Datta, kappa = 1/M_Pl^2
    G_strong = c / Lam2                               # independent strong coupling
    return {
        "G_crit": G_crit,
        "ratio_grav": sp.simplify(G_grav / G_crit),   # ~ (Lambda/M_Pl)^2
        "ratio_strong": sp.simplify(G_strong / G_crit),  # ~ c, O(1)
        "c_for_criticality": sp.simplify(sp.solve(sp.Eq(G_strong / G_crit, 1), c)[0]),
    }


def mvac_r0(g):
    """m_vac*r0 in the framework's substrate relations: M = g v, sigma = 2 pi v^2,
    r0 = 1.1/sqrt(sigma) (the QCD Sommer value r0*sqrt(sigma) ~= 1.1). Then m_vac*r0 = 1.1 g/sqrt(2 pi)."""
    return 1.1 * g / sp.sqrt(2 * sp.pi)


def report():
    import sympy as sp
    G_crit, I0_0 = njl_critical_coupling()
    print("NJL critical coupling (gap turns on iff G_S > G_crit):")
    print(f"    I0(M->0) = {I0_0}")
    print(f"    G_crit   = {G_crit}\n")

    fk = fork()
    print("(A) gravitational Hehl-Datta  G_S = 3 kappa/16 = 3/(16 M_Pl^2):")
    print(f"    G_S / G_crit = {fk['ratio_grav']}   (~ (Lambda/M_Pl)^2)")
    print("(B) independent strong        G_S = c/Lambda^2:")
    print(f"    G_S / G_crit = {fk['ratio_strong']}   -> criticality at c = {fk['c_for_criticality']}\n")

    # --- numbers ---
    Lam_GeV, MPl_GeV, Nc = 1.0, 2.435e18, 3            # Lambda ~ 1 GeV, reduced Planck mass
    ratio_grav = float(fk["ratio_grav"].subs({sp.Symbol("Lambda", positive=True): Lam_GeV,
                                              sp.Symbol("M_Pl", positive=True): MPl_GeV,
                                              sp.Symbol("N_c", positive=True): Nc}))
    c_crit = float(fk["c_for_criticality"].subs(sp.Symbol("N_c", positive=True), Nc))
    print(f"numbers (Lambda={Lam_GeV} GeV, M_Pl={MPl_GeV:.3g} GeV, N_c={Nc}):")
    print(f"    (A) G_S/G_crit = {ratio_grav:.2e}   -> sub-critical by ~37 orders -> M = 0 EXACTLY")
    print(f"    (B) criticality at c = {c_crit:.2f}; gap turns on for c > that\n")

    g_sym = sp.symbols("g", positive=True)
    g_for_3 = float(sp.solve(sp.Eq(mvac_r0(g_sym), 10), g_sym)[0])
    g_qcd = 3.76                                        # M/f_pi ~ 350/93
    print("coupling strength implied (m_vac*r0 = 1.1 g/sqrt(2 pi), g = M/f_pi):")
    print(f"    QCD-like g = {g_qcd:.1f}  -> m_vac*r0 = {float(mvac_r0(g_qcd)):.2f}  (~1 bound level)")
    print(f"    three generations (m_vac*r0 = 10) need g = M/f_pi ~ {g_for_3:.0f}  "
          f"(~{g_for_3/g_qcd:.0f}x deeper into the broken phase than QCD)\n")

    print("VERDICT: gravitational torsion gives the gap NO chance (sub-critical by ~1e-37): no mass,")
    print("no torsiton, no IR matter. The framework's matter REQUIRES torsion to be an independent")
    print("strong field, G_S ~ c/Lambda^2 with c well above critical -- decoupled from Newton's G.")
    print("That is a dynamical postulate ('torsion is its own thing'), not a result of Einstein-Cartan.")


if __name__ == "__main__":
    report()
