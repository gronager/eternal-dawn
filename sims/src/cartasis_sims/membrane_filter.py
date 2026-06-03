r"""Tier 3: the membrane filter and the dark-to-baryon ratio f ~ 1/6.

If dark matter is the parent matter the bounce membrane *rejects* -- felt
gravitationally but never crossing into our plasma (Ch.5, Ch.8 D) -- then the
dark-to-baryon ratio fixes the membrane PASS-FRACTION

    f = omega_b / (omega_b + omega_c),     omega_c/omega_b ~ 5.4  ->  f ~ 0.157 ~ 1/6.

This module asks what Einstein--Cartan microphysics sets f. The honest answer is a
NARROWING, not a parameter-free derivation, and it is stated as such.

THE MODEL. The membrane is a torsion barrier: the spin--spin (Hehl--Datta)
repulsion that halts the collapse is also a potential barrier Delta that infalling
matter must overcome to cross. For a thermal gas at the bounce temperature T the
pass-fraction is then a barrier penetration / occupation factor controlled by the
single dimensionless ratio x = Delta/T:

    Boltzmann (Maxwell tail above a barrier):  f = exp(-x),
    Fermi--Dirac (states above a gap)       :  f = 1 / (exp(x) + 1).

Inverting (verified with SymPy in the tests):
    f = 1/6  <=>  x = ln 6 = 1.79  (Boltzmann)   or   x = ln 5 = 1.61  (Fermi),
and the measured f = 0.157 sits at x ~ 1.68--1.85.

THE ANCHOR (why this is not numerology). The membrane is, by definition (Ch.2), the
surface where the torsion energy density equals the matter energy density --
repulsion balances gravity. Per particle that means the barrier and the thermal
scale are COMPARABLE there, x = Delta/T ~ O(1): not 10^-9, not 10^9. An O(1) barrier
gives f in the 0.08--0.37 decade (x in [1, 2.5]); the observed 1/6 lands squarely
inside, at x ~ 1.8 -- the barrier ~80% above the mean thermal energy, i.e. torsion
*just* dominating, exactly the membrane condition.

WHAT THIS BUYS, AND WHAT IT DOES NOT. It converts f from a free fit in [0,1] into an
O(1) barrier ratio the bounce supplies automatically -- so f ~ 1/6 is *natural*, not
tuned. It does NOT yet derive 0.157 to three digits: pinning x = 1.8 rather than 1.0
or 2.5 needs the full inhomogeneous bounce profile (a Tier-3 computation). The
bounce temperature scale here is ~10^7 GeV (cve_filter.py).
"""

from __future__ import annotations

import math

from . import constants as k

# bounce-energy scale (GeV); the thermal scale T in x = Delta/T
BOUNCE_T_GEV = 1.05e7

# an "O(1) barrier" range for the torsion-balance membrane
X_NATURAL_LO = 1.0
X_NATURAL_HI = 2.5


def pass_fraction_boltzmann(x: float) -> float:
    """Maxwell-tail pass-fraction over a barrier x = Delta/T: f = exp(-x)."""
    return math.exp(-x)


def pass_fraction_fermi(x: float) -> float:
    """Fermi--Dirac occupation above a gap x = Delta/T: f = 1/(exp(x)+1)."""
    return 1.0 / (math.exp(x) + 1.0)


def barrier_ratio_for_f(f: float, kind: str = "boltzmann") -> float:
    """Invert the pass-fraction: the barrier ratio x = Delta/T giving pass-fraction f.
    Boltzmann: x = -ln f ; Fermi: x = ln((1-f)/f)."""
    if not 0.0 < f < 1.0:
        raise ValueError("f must be in (0,1)")
    if kind == "boltzmann":
        return -math.log(f)
    if kind == "fermi":
        return math.log((1.0 - f) / f)
    raise ValueError("kind must be 'boltzmann' or 'fermi'")


def observed_f(omega_b: float = k.Omega_b_hsq,
               omega_c: float = k.Omega_c_hsq) -> float:
    """The membrane pass-fraction read from the Planck densities: f = ob/(ob+oc)."""
    return omega_b / (omega_b + omega_c)


def natural_band(kind: str = "boltzmann",
                 x_lo: float = X_NATURAL_LO, x_hi: float = X_NATURAL_HI):
    """Range of f produced by an O(1) torsion barrier x in [x_lo, x_hi].
    Returns (f_at_x_hi, f_at_x_lo) = (f_min, f_max) since f decreases in x."""
    fn = pass_fraction_boltzmann if kind == "boltzmann" else pass_fraction_fermi
    return fn(x_hi), fn(x_lo)


def is_f_natural(f: float = None, kind: str = "boltzmann") -> bool:
    """Does the observed f fall inside the O(1)-barrier band (i.e. is it natural,
    not fine-tuned to ~0 or ~1)?"""
    f = observed_f() if f is None else f
    lo, hi = natural_band(kind)
    return lo <= f <= hi
