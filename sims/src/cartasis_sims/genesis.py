r"""OGU genesis: the rotating-bounce vortex, what it entrains, and the OGU mass.

Does the spinning Kerr--Cartan bounce PUMP void matter inward, and what patch mass
does it entrain before freezing? Two pieces.

1. The pump (spin window). A rotating bounce drags frames (Lense--Thirring,
   Omega_drag ~ 2GJ/c^2 r^3): the throat is a vortex that spins up and draws in the
   surrounding void, the rotational analogue of Bondi capture. The spin must clear
   TWO gates. Below f_min = eta_min/(C * (Omega_bounce/T)) the chiral-vortical
   asymmetry is too small to survive annihilation -- a sterile seed that pumps
   nothing (Ch.5 spin filter). Near extremal (f -> 1) the centrifugal barrier
   chokes radial inflow (matter circularizes instead of falling in). So the
   entrainment efficiency is peaked at intermediate spin -- a genesis sweet spot,
   and the "productivity" factor of the observed-spin distribution.

2. The mass (horizon-limited growth). Once pumping, the seed grows by runaway
   accretion (extruder/growth.py) and entrains the void within its causal reach,
   R ~ c t. A patch of radius ct holds at most the mass whose Schwarzschild radius
   equals ct,

       M_horizon(t) = c^3 t / (2 G),

   so if the void is at least critically dense (rho_void >= rho_crit(t) =
   3/(8 pi G t^2) -- a "condensed" void), the patch fills and COLLAPSES to that
   horizon mass: the OGU becomes a horizon-scale hole and

       M_OGU ~ c^3 t / (2 G),

   set by its AGE, not the void density. (For a sub-critical void the OGU is
   density-limited, M ~ rho_void (4/3) pi (ct)^3, and smaller.) So the OGU mass
   measures the supraverse age: M_OGU ~ our mass at one Hubble time, ~1e65 kg at
   ~1e12 Hubble times. Its sibling count is N_sib ~ epsilon M_OGU/M_vis, so a
   1e65 kg OGU has ~1e11 universes like ours -- yet we are still BHU1
   (population.py: the depth is set by descendants pinned at M_vis, not the OGU).

The headline relation M_OGU ~ c^3 t/(2G) is robust (it is just the horizon-mass /
"universe = its own Schwarzschild radius" identity); the spin sweet-spot shape is
an illustrative parametrization of the frame-dragging-vs-choking competition.
"""

from __future__ import annotations

import math

from . import constants as k
from . import cve_filter as cve

FOUR_THIRDS_PI = 4.0 / 3.0 * math.pi
HUBBLE_TIME_S = 4.35e17          # ~ our age, an illustrative supraverse-age unit
M_OUR = 9.2e52                   # our universe's mass [kg] ~ M_vis


# --- the entrained mass: horizon-limited growth ----------------------------

def horizon_mass(t_s: float) -> float:
    """Mass whose Schwarzschild radius equals the causal radius c t: c^3 t/(2G)."""
    return k.c**3 * t_s / (2.0 * k.G)


def void_critical_density(t_s: float) -> float:
    """Critical density at age t (H ~ 1/t): rho_crit = 3/(8 pi G t^2)."""
    return 3.0 / (8.0 * math.pi * k.G * t_s**2)


def density_limited_mass(t_s: float, rho_void: float) -> float:
    """Void matter within the causal patch: rho_void (4/3) pi (c t)^3."""
    return rho_void * FOUR_THIRDS_PI * (k.c * t_s)**3


def ogu_mass(t_s: float, rho_void: float | None = None) -> float:
    """OGU mass at age t. Horizon-limited (=c^3 t/2G) if the void is at least
    critically dense; otherwise density-limited. rho_void=None -> condensed
    (horizon-limited)."""
    Mh = horizon_mass(t_s)
    if rho_void is None:
        return Mh
    return min(density_limited_mass(t_s, rho_void), Mh)


def age_for_mass(M_ogu: float) -> float:
    """Supraverse age implied by a horizon-mass OGU: t = 2 G M/c^3."""
    return 2.0 * k.G * M_ogu / k.c**3


def n_siblings(M_ogu: float, M_vis: float = M_OUR, epsilon: float = 0.1,
               f_clean: float = 1.0) -> float:
    """Viable BHU1 children an OGU spawns: ~ epsilon f_clean M_ogu/M_vis."""
    return epsilon * f_clean * M_ogu / M_vis


# --- the pump: the spin window and entrainment sweet spot ------------------

def spin_floor(eta_min: float = 2.0e-11, C: float = 10.0) -> float:
    """f_min = eta_min/(C * Omega_bounce/T): below this the seed is sterile.
    f is the spin fraction omega/Omega_bounce (<=1, sub-extremal)."""
    return eta_min / (C * cve.vorticity_over_T(spin_fraction=1.0))


def eta_from_spin(f: float, C: float = 10.0) -> float:
    """Chiral-vortical asymmetry of a seed at spin fraction f: eta = C f (Omega/T)."""
    return C * f * cve.vorticity_over_T(spin_fraction=1.0)


def entrainment_efficiency(f: float, eta_min: float = 2.0e-11,
                           C: float = 10.0) -> float:
    """Illustrative vortex pump efficiency vs spin fraction f in [0,1]:
    frame-dragging rises (~f^2), centrifugal choking cuts in near extremal
    (~1-f^2), and survival gates it below f_min. Peaks at the genesis sweet spot."""
    fmin = spin_floor(eta_min, C)
    if f <= fmin or f >= 1.0:
        return 0.0
    return f**2 * (1.0 - f**2)
