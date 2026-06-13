r"""The sign: does the curvature of a bounce melt the chiral condensate, or catalyse it?

In curved space the squared Dirac operator is (i D)^2 = -nabla^2 + R/4 + M^2 (Lichnerowicz), so the
fermion loop sees an EFFECTIVE mass M_eff^2 = M^2 + R/4. The NJL gap equation becomes

        1 = 8 G N_c I0(M^2 + R/4),     I0(x) = (1/16 pi^2)[Lambda^2 - x ln(1 + Lambda^2/x)].

Increasing R raises M_eff^2 and SUPPRESSES the loop -> the condensate is restored (mass melts).
The condensate survives only up to a critical curvature; equating the M=0 onset to the flat-space gap
(1 = 8 G N_c I0(M_flat^2)) gives the clean result

        R_crit = 4 M_flat^2     (curvature radius L_crit = 1/sqrt(R_crit) = 1/(2 M_flat) ~ Compton/2).

So: R>0 RESTORES (the bounce melts the mass -- the good sign), R<0 catalyses; and restoration sets in
when the curvature radius reaches the constituent Compton wavelength -- the SAME scale as the Unruh
melt and the holographic void seed. Low curvature (flat, expanded) -> condensate on -> fermions exist.
This module shows the sign and R_crit symbolically and checks the bounce and a neutron star.
"""
from __future__ import annotations

import sympy as sp


def I0(x, Lam2):
    """NJL loop int d^4k_E/(2pi)^4 1/(k^2+x), sharp cutoff Lambda^2, as a function of x = M_eff^2."""
    return (Lam2 - x * sp.log(1 + Lam2 / x)) / (16 * sp.pi**2)


def critical_curvature():
    """R_crit where the curved gap loses its broken solution. Returns the symbolic identity."""
    Mf, R = sp.symbols("M_flat R", positive=True)
    # flat gap fixes 8 G N_c I0(M_flat^2) = 1; restoration onset is M=0 with M_eff^2 = R/4:
    # 8 G N_c I0(R/4) = 1 = 8 G N_c I0(M_flat^2)  ->  R/4 = M_flat^2.
    return sp.Eq(R, 4 * Mf**2)


def gap_mass(R_val, Mflat=0.35):
    """Curved-gap dynamical mass at Ricci scalar R (GeV^2). Because the gap 8 G Nc I0(M^2+R/4)=1 with
    the SAME coupling that fixes the flat gap 8 G Nc I0(M_flat^2)=1, and I0 is strictly monotonic,
    the loop arguments must match: M^2 + R/4 = M_flat^2. So EXACTLY

        M(R) = sqrt(M_flat^2 - R/4)   for R < 4 M_flat^2,   else 0 (restored).

    A mean-field second-order melting: M shrinks as the curvature grows, vanishing at R_crit=4 M_flat^2."""
    import math
    m2 = Mflat**2 - R_val / 4.0
    return math.sqrt(m2) if m2 > 0 else 0.0


def report():
    print("Sign of the curvature effect on the chiral condensate (Lichnerowicz M_eff^2 = M^2 + R/4):\n")
    print(f"  symbolic onset:  {critical_curvature()}   ->  R_crit = 4 M_flat^2")
    print("  => R>0 RESTORES (melts the mass), R<0 catalyses. The bounce has R>0: the GOOD sign.\n")

    Mflat = 0.35           # GeV, a QCD-like constituent/dynamical mass
    print(f"  curved gap, M_flat = {Mflat} GeV  (R_crit = {4*Mflat**2:.3f} GeV^2):")
    print("    R (GeV^2)   M_dyn (GeV)   phase")
    for R in (0.0, 0.1, 0.3, 0.49, 0.6, 1.0):
        M = gap_mass(R, Mflat=Mflat)
        phase = "broken (mass ON)" if M > 1e-6 else "RESTORED (mass off)"
        print(f"      {R:>5.2f}      {M:>8.4f}     {phase}")

    # the restoration scale as a length, vs the Unruh/holographic scale and a neutron star
    GeV_per_invfm = 0.1973                                  # hbar c
    L_crit_fm = GeV_per_invfm / (2 * Mflat)                 # 1/(2 M_flat) in fm
    print(f"\n  restoration curvature radius L_crit = 1/(2 M_flat) = {L_crit_fm:.3f} fm "
          f"(~ constituent Compton)")
    print("    -> same scale as the Unruh melt (~0.03-0.16 fm) and the holographic void seed (~0.4 fm).")
    # neutron star: R ~ 8 pi G rho /c^2 ; curvature radius vs L_crit
    G, c = 6.674e-11, 2.998e8
    rho_ns = 1.0e18
    R_ns = 8 * 3.141592653589793 * G * rho_ns / c**2        # 1/m^2
    L_ns_fm = (1 / R_ns**0.5) * 1e15
    print(f"  neutron star (rho~1e18): curvature radius ~ {L_ns_fm:.1e} fm  >>  L_crit "
          f"-> mass safe, no melt, no bounce.")
    print("\n  VERDICT: the bounce's positive curvature MELTS the condensate (fermions cannot crystallise")
    print("  at the bounce); flat/expanded space (R << M^2) keeps it broken (fermions exist). The sign")
    print("  is the one the framework needs -- a real, grown, low-curvature world, not a Boltzmann brain.")


if __name__ == "__main__":
    report()
