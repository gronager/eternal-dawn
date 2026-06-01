r"""The dark sector as look-back into the parent: gravity felt before light seen.

This ties together the three things the framework says about the parent's ongoing
accretion through the Cartasis membrane (Ch.6), and gets the full dynamic right.

The membrane is our PAST cosmological horizon -- our Big-Bang surface (Ch.6: CMB =
parent Hawking radiation at that surface). Matter the parent injects through it
therefore arrives in our DEEP PAST, not our present: from inside, the parent's
"now" is our "earliest visible". Three consequences, all computed from flat LambdaCDM:

1. We see further back every year. The comoving PARTICLE horizon grows monotonically:
   chi_p(t) = c \int_0^t dt'/a(t'). Today chi_p ~ 47 Gly; it rises toward an asymptote
   chi_p(inf) ~ 64 Gly as the de Sitter future runs. So if today we look back ~13.8 Gyr
   (z~1100 at the membrane), in +1 Gyr we look back FURTHER -- the horizon grows ~1 Gly.
   New parent-injected matter, arriving at the receding past horizon, becomes visible
   only as the horizon expands out to it. (We will ever see only ~74% of the way to the
   asymptote -- the rest is still climbing toward us.)

2. But we FEEL its gravity before we SEE it (dark matter). The parent's mass sits at /
   beyond the membrane and acts gravitationally on our interior NOW -- it is already in
   our metric (the projected mass curves our space), while its LIGHT is still climbing
   in from the past horizon. So the dark-matter gravitational signature LEADS the
   electromagnetic visibility: we measure the pull of parent matter we cannot yet see.
   Dark matter is exactly "gravity arriving ahead of light" from the membrane.

3. The dark sector is a RATE, not a relic. Dark matter ~ the projected parent mass
   (its accretion history); dark energy ~ the time-derivative of the parent's growth
   (the membrane geometry changing). Both evolve as the parent's meal evolves -- the
   handle on the sigma_8 / S_8 and DESI w0-wa hints (Ch.6).

Units: flat LambdaCDM, Planck 2018; horizons returned in giga-light-years.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from . import constants as k

GLY = 9.4607304725808e24 * 1.0e9 / 1.0e9    # 1 Gly in metres (= c * Gyr)
GLY_M = k.c * (k.year * 1.0e9)              # metres in a giga-light-year
GYR_S = k.year * 1.0e9                      # seconds in a Gyr


def _H(a, Om=k.Omega_m, OL=k.Omega_Lambda):
    """Flat LambdaCDM expansion rate H(a) [s^-1]."""
    return k.H0 * np.sqrt(Om / a**3 + OL)


def cosmic_age(a=1.0, Om=k.Omega_m, OL=k.Omega_Lambda):
    """Cosmic time at scale factor a [s]."""
    return quad(lambda ap: 1.0 / (ap * _H(ap, Om, OL)), 1e-8, a, limit=200)[0]


def particle_horizon(a=1.0, Om=k.Omega_m, OL=k.Omega_Lambda):
    """Comoving particle horizon chi_p(a) = c int_0^a da'/(a'^2 H) [m]: the farthest
    BACK we can see -- the radius of the visible interior, hence of the parent meal
    we have received light from."""
    return k.c * quad(lambda ap: 1.0 / (ap**2 * _H(ap, Om, OL)), 1e-7, a, limit=200)[0]


def event_horizon(a=1.0, Om=k.Omega_m, OL=k.Omega_Lambda):
    """Comoving event horizon [m]: parent matter injected beyond this comoving radius
    will NEVER become visible (its light never reaches us)."""
    return k.c * quad(lambda ap: 1.0 / (ap**2 * _H(ap, Om, OL)), a, 1e3, limit=200)[0]


def scale_factor_at_age(t_s, Om=k.Omega_m, OL=k.Omega_Lambda):
    """Invert cosmic_age: the scale factor at cosmic time t [s]."""
    from scipy.optimize import brentq
    return brentq(lambda a: cosmic_age(a, Om, OL) - t_s, 1e-6, 1e3)


def lookback_gain(delta_gyr=1.0, Om=k.Omega_m, OL=k.Omega_Lambda):
    """How much FURTHER back (in Gly of comoving particle horizon) we can see after
    waiting delta_gyr giga-years. Positive: the past horizon recedes, revealing more
    of the parent's injected matter."""
    t0 = cosmic_age(1.0, Om, OL)
    a1 = scale_factor_at_age(t0 + delta_gyr * GYR_S, Om, OL)
    return (particle_horizon(a1, Om, OL) - particle_horizon(1.0, Om, OL)) / GLY_M


def horizon_gly(a=1.0):
    """Particle horizon in giga-light-years (today ~47)."""
    return particle_horizon(a) / GLY_M


def asymptotic_horizon_gly():
    """The comoving particle horizon as a -> inf (de Sitter future), in Gly (~64):
    the most of the parent's history we will EVER see."""
    return particle_horizon(1e3) / GLY_M


def fraction_seen_now():
    """Fraction of the ever-visible parent history we have received light from today."""
    return horizon_gly(1.0) / asymptotic_horizon_gly()
