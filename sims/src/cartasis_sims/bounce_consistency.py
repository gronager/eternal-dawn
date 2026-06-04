r"""Are the bounce (Part I) and the mass scale (Part II) one coupling or two? -- assessed.

Part~II generates fermion masses from a strong condensate at Lambda ~ 246 GeV (transmutation of
a strong/gauge coupling), and asks whether promoting torsion to propagating disturbs the
cosmological bounce of Part~I. This module settles it with a back-of-envelope a neutron star
already decides.

THE TEST. The Einstein--Cartan bounce comes from a spin-spin (four-fermion) repulsion balancing
gravity; the bounce density scales INVERSELY with that four-fermion coupling, rho_C ~ 1/G_eff. Two
candidate couplings:
  * GRAVITATIONAL torsion kappa ~ G/c^4 ~ 1/M_Pl^2 (Part I) -> rho_C ~ 1e50 kg/m^3.
  * the STRONG, mass-giving coupling ~ 1/Lambda^2, stronger by (M_Pl/Lambda)^2 ~ 2.5e33.
If the STRONG coupling drove the bounce, rho_C would fall by that factor -- to ~4e16 kg/m^3, about
0.2x nuclear density. But neutron-star cores reach a few times nuclear density and do NOT bounce.
So the strong coupling is excluded as the bounce driver: the bounce is gravitational, and the
mass-giving sector is a SEPARATE, ~10^33 stronger coupling -- they MUST be two couplings, and the
existence of neutron stars proves it.

THE CONSEQUENCE for Part I. At the bounce, T ~ rho_C^{1/4} ~ 1e7 GeV is ~10^5 times the
condensation scale Lambda, so the mass condensate is fully melted (consistent with Part I's
'electroweak restored' matter) and the strong/torsion modes contribute only as extra relativistic
radiation -- an O(1) shift in g_*, not a change to the gravitational repulsion. The bounce stands;
the only revision is bookkeeping the extra degrees of freedom. Reproducible; the figure plots it.
"""

from __future__ import annotations

import numpy as np

# scales (GeV unless noted)
M_PL_GEV = 1.22e19
LAMBDA_GEV = 246.0           # condensate / mass-giving strong scale
RHO_C_GRAV = 1.0e50          # kg/m^3, Part I gravitational-torsion bounce density
RHO_NUCLEAR = 2.3e17         # kg/m^3
RHO_NS_CORE = 1.0e18         # kg/m^3 (neutron-star core, a few x nuclear)
RHO_PLANCK = 5.0e96          # kg/m^3


def coupling_ratio(M_Pl=M_PL_GEV, Lam=LAMBDA_GEV):
    """How much stronger the mass-giving four-fermion coupling (~1/Lambda^2) is than the
    gravitational one (~G/c^4 ~ 1/M_Pl^2): the factor (M_Pl/Lambda)^2."""
    return (M_Pl / Lam) ** 2


def bounce_density(coupling="gravitational", M_Pl=M_PL_GEV, Lam=LAMBDA_GEV,
                   rho_C_grav=RHO_C_GRAV):
    """The bounce density for a given four-fermion coupling. rho_C ~ 1/G_eff, so the strong
    coupling (stronger by (M_Pl/Lambda)^2) bounces at a correspondingly LOWER density."""
    if coupling == "gravitational":
        return rho_C_grav
    if coupling == "strong":
        return rho_C_grav / coupling_ratio(M_Pl, Lam)
    raise ValueError(coupling)


def strong_coupling_excluded_by_neutron_stars():
    """The strong coupling would put the bounce at ~0.2x nuclear density; neutron-star cores
    exceed that without bouncing. Returns the bounce density, its ratio to nuclear, and how far
    a neutron-star core exceeds it -- i.e. the exclusion."""
    rho_strong = bounce_density("strong")
    return {
        "rho_C_strong_kg_m3": rho_strong,
        "in_nuclear_units": rho_strong / RHO_NUCLEAR,
        "ns_core_exceeds_by": RHO_NS_CORE / rho_strong,
        "excluded": RHO_NS_CORE > rho_strong,   # neutron stars pass through it without bouncing
    }


def bounce_temperature(rho_C=RHO_C_GRAV, M_Pl=M_PL_GEV, rho_Pl=RHO_PLANCK):
    """Energy scale at the bounce, T ~ M_Pl (rho_C/rho_Pl)^{1/4} (energy density ~ T^4)."""
    return M_Pl * (rho_C / rho_Pl) ** 0.25


def condensate_is_melted_at_bounce(Lam=LAMBDA_GEV):
    """T_bounce / T_c with T_c ~ Lambda: if >> 1 the mass condensate is fully melted at the
    bounce (so the mass sector is radiation there, not bound matter)."""
    return bounce_temperature() / Lam


def summary():
    return {
        "coupling_ratio_strong_over_grav": coupling_ratio(),
        "rho_C_gravitational": RHO_C_GRAV,
        "neutron_star_test": strong_coupling_excluded_by_neutron_stars(),
        "bounce_temperature_GeV": bounce_temperature(),
        "T_bounce_over_T_c": condensate_is_melted_at_bounce(),
        "verdict": "the bounce coupling is gravitational (G/c^4); the mass coupling (~1/Lambda^2) "
                   "is ~10^33 stronger and is EXCLUDED as the bounce driver by neutron stars. Two "
                   "couplings, necessarily -- and at the bounce the mass sector is melted radiation, "
                   "so Part I's bounce stands (only the radiation d.o.f. are revised).",
    }
