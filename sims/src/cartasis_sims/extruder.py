r"""Tier 1c: the extruder -- the entropy budget of the Einstein--Cartan bounce.

    "Universes are not born, they are extruded."

The bounce does not *sort* a pre-existing population (a filter); it compresses
parent material to the Cartasis density rho_C and extrudes a new expanding
interior through it, like material forced through a die. Sphalerons conserve
B-L, so the parent's net B-L passes through as a protected remnant and the baby
inherits

    eta_baby = eta_parent / D,

where D = S_after / S_before >= 1 is the entropy multiplied across the crossing
(`cve_filter.inherited_eta`). The whole magnitude-and-lineage question reduces to
this one number D. This module computes it from the bounce dynamics rather than
guessing it.

Two anchors bracket D (see `cve_filter.dilution_horizon` for the upper one):

  * D = 1            -- an adiabatic (reversible) extrusion conserves entropy and
                        passes eta through exactly: pure inheritance.
  * D = S_BH/S_rad   -- the interior thermalizes to the parent horizon entropy
        ~ 1e49          (M ~ 1e53 kg): eta crushed to ~1e-59, which we DO NOT see.

So which is it? Entropy is generated only by *irreversibility*. For a perfect
fluid the homogeneous bounce is exactly reversible (D = 1). Switch on dissipation
as a bulk viscosity zeta = zeta_tilde * s (zeta_tilde dimensionless, s the
entropy density). The Eckart entropy production of a bulk-viscous FRW fluid is

    d(ln S)/dtau = zeta (3H)^2 / (s T) = 9 zeta_tilde H^2 / T          (natural units),

so, integrating across the bounce,

    ln D = 9 zeta_tilde * (Omega_bounce / T_bounce) * J,
    J = \int H_dim^2 a  dtau_dim   (a pure number from the dimensionless bounce),

with Omega_bounce = sqrt(8 pi G rho_C / 3) and T_bounce the plasma temperature.
The controlling prefactor Omega_bounce/T_bounce is the SAME ~1e-11 small
parameter that limits the chiral-vortical estimate (`cve_filter`): the bounce is
fast compared to the thermal time, so it has no time to generate entropy. The
result is that the extruder is adiabatic to ~1 part in 1e11 even for O(1)
viscosity -- and, since a radiation (conformal) fluid has zeta = 0 identically,
to begin with there is no bulk-viscous entropy at all. Inheritance is robust;
the horizon-saturating limit is dynamically out of reach for a fast bounce.

Caveats (uncomputed, flagged): this is the homogeneous channel. Inhomogeneous
shear, shocks, and gravitational particle production at the bounce are separate
entropy sources not captured by an Eckart bulk viscosity; the timescale argument
suppresses anything sourced by gradients x (bounce time) by the same ~1e-11
unless the gradients are trans-thermal. Tier 1 (inhomogeneous collapse) is where
those are tested.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import bounce as bnc
from . import cve_filter as cve


@dataclass
class ExtruderResult:
    w: float
    rho_C: float
    J: float                 # dimensionless entropy-production integral
    adiabaticity: float      # Omega_bounce / T_bounce  (~1e-11)
    T_bounce_GeV: float
    Omega_bounce_GeV: float


def temperature_over_Tbounce(a: np.ndarray) -> np.ndarray:
    """Radiation plasma temperature in units of T_bounce: T/T_b = a_min/a = 1/a
    (a is in units of a_min = 1; for a conformal fluid s ~ T^3 ~ a^-3)."""
    return 1.0 / a


def entropy_production_integral(w: float = 1.0 / 3.0, a_init: float = 6.0,
                                t_max: float = 30.0,
                                n_points: int = 8000) -> float:
    r"""J = \int H_dim^2 a dtau_dim across the bounce (dimensionless).

    H_dim and a come straight from the dimensionless modified-Friedmann bounce
    (`bounce.simulate_bounce`); the integrand H^2 a is the dimensionless kernel
    of the Eckart entropy production once temperature (T ~ 1/a) is folded in.
    """
    sol = bnc.simulate_bounce(w=w, a_init=a_init, t_max=t_max, n_points=n_points)
    integrand = sol.H**2 * sol.a
    return float(np.trapezoid(integrand, sol.t))


def analyze(w: float = 1.0 / 3.0, rho_C: float = 1.0e50,
            a_init: float = 6.0, t_max: float = 30.0,
            n_points: int = 8000) -> ExtruderResult:
    """Bundle the dimensionless integral with the physical bounce scales."""
    J = entropy_production_integral(w, a_init, t_max, n_points)
    T_b = cve.bounce_temperature_GeV(rho_C)
    Om = cve.bounce_rate_GeV(rho_C)
    return ExtruderResult(w=w, rho_C=rho_C, J=J, adiabaticity=Om / T_b,
                          T_bounce_GeV=T_b, Omega_bounce_GeV=Om)


def ln_dilution(zeta_tilde: float, res: ExtruderResult | None = None,
                **kw) -> float:
    """ln D = 9 zeta_tilde (Omega_bounce/T_bounce) J for a bulk viscosity
    zeta = zeta_tilde * s. Conformal (radiation) matter has zeta_tilde = 0."""
    if res is None:
        res = analyze(**kw)
    return 9.0 * zeta_tilde * res.adiabaticity * res.J


def dilution_factor(zeta_tilde: float, res: ExtruderResult | None = None,
                    **kw) -> float:
    """D = exp(ln D) >= 1, the entropy multiplied across the extrusion."""
    return float(np.exp(ln_dilution(zeta_tilde, res, **kw)))


def inherited_eta(zeta_tilde: float, eta_parent: float = cve.ETA_OBS,
                  res: ExtruderResult | None = None, **kw) -> float:
    """eta_baby = eta_parent / D for the dynamically-computed dilution D."""
    return cve.inherited_eta(eta_parent, dilution_factor(zeta_tilde, res, **kw))
