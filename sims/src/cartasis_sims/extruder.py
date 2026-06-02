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


# ----------------------------------------------------------------------------
# Inhomogeneous shear channel.
#
# The homogeneous channel above used a bulk viscosity (sourced by the expansion
# theta = 3H). The worry is the *shear*: an anisotropic/inhomogeneous collapse
# carries a shear scalar sigma whose energy density redshifts as a stiff fluid,
# rho_sigma ~ a^-6, so it grows FASTER than radiation (a^-4) and comes to
# dominate the approach to the bounce (the BKL/Mixmaster regime). In plain GR
# that drives sigma^2/T -> infinity and the shear-viscous entropy production
# DIVERGES at the singularity. The question is what the Einstein--Cartan bounce
# does to it.
#
# Minimal model: a Bianchi-I-like mean scale factor with the torsion cap on the
# TOTAL effective density (dimensionless rho_C = 1, 8 pi G/3 = 1, a_min = 1):
#
#     rho_tot(a) = (1 - f_sigma) a^-n + f_sigma a^-6,    n = 3(1+w),
#     H^2 = rho_tot (1 - rho_tot),
#
# where f_sigma is the shear fraction of the density AT the bounce (a = 1):
# f_sigma = 0 is the isotropic bounce, f_sigma -> 1 a maximally shear-dominated
# one. The shear-viscous (Eckart) entropy production is d ln S/dtau = 2 eta~ sigma^2/T
# with sigma^2 proportional to rho_sigma = f_sigma a^-6, so
#
#     ln D_shear = 2 eta~ (Omega_bounce/T_bounce) J_shear,
#     J_shear = int rho_sigma(a) a dtau   (a pure number, O(1-10)).
#
# The SAME prefactor Omega/T ~ 1e-11 appears, because the torsion cap bounds the
# shear at the bounce scale (sigma <~ Omega): the bounce is fast in thermal
# units whatever drives it. So even a maximally shear-dominated, KSS-viscous
# bounce stays adiabatic, and the same torsion that removes the singularity also
# tames the BKL entropy divergence.
# ----------------------------------------------------------------------------

def _rho_tot(a, f_sigma: float, n: float):
    return (1.0 - f_sigma) * a**(-n) + f_sigma * a**(-6.0)


def _rho_tot_prime(a, f_sigma: float, n: float):
    return (-n * (1.0 - f_sigma) * a**(-n - 1.0)
            - 6.0 * f_sigma * a**(-7.0))


def simulate_anisotropic_bounce(f_sigma: float = 0.5, w: float = 1.0 / 3.0,
                                tau_max: float = 30.0, n_points: int = 8000):
    """Integrate the shear-augmented torsion bounce outward from a_min = 1.

    Returns (tau, a, H, rho_sigma) on the expansion branch; the collapse branch
    is its mirror image, so entropy integrals over the full bounce are 2x this.
    """
    from scipy.integrate import solve_ivp

    n = 3.0 * (1.0 + w)

    def rhs(tau, y):
        a, v = y
        rt = _rho_tot(a, f_sigma, n)
        F = rt * (1.0 - rt)
        Fp = _rho_tot_prime(a, f_sigma, n) * (1.0 - 2.0 * rt)
        a_dd = a * F + 0.5 * a**2 * Fp
        return [v, a_dd]

    sol = solve_ivp(rhs, (0.0, tau_max), [1.0, 0.0],
                    t_eval=np.linspace(0.0, tau_max, n_points),
                    rtol=1e-10, atol=1e-12, method="DOP853")
    a = sol.y[0]
    H = sol.y[1] / a
    rho_sigma = f_sigma * a**(-6.0)
    return sol.t, a, H, rho_sigma


def shear_entropy_integral(f_sigma: float = 0.5, w: float = 1.0 / 3.0,
                           tau_max: float = 30.0, n_points: int = 8000) -> float:
    r"""J_shear = \int rho_sigma a dtau over the full bounce (= 2x expansion)."""
    tau, a, H, rho_sigma = simulate_anisotropic_bounce(f_sigma, w, tau_max,
                                                        n_points)
    return 2.0 * float(np.trapezoid(rho_sigma * a, tau))


def shear_ln_dilution(eta_tilde: float, f_sigma: float = 0.5,
                      rho_C: float = 1.0e50, w: float = 1.0 / 3.0,
                      **kw) -> float:
    """ln D from shear viscosity eta_s = eta_tilde * s across the bounce."""
    T_b = cve.bounce_temperature_GeV(rho_C)
    Om = cve.bounce_rate_GeV(rho_C)
    J = shear_entropy_integral(f_sigma, w, **kw)
    return 2.0 * eta_tilde * (Om / T_b) * J


def shear_dilution_factor(eta_tilde: float, f_sigma: float = 0.5, **kw) -> float:
    return float(np.exp(shear_ln_dilution(eta_tilde, f_sigma, **kw)))

