r"""Tying the spin story together, and anchoring eta -- what is true and observable.

Two quantities that the framework keeps strictly SEPARATE are easy to conflate:

  CHIRALITY (the sign of eta: matter vs antimatter) is INHERITED down a lineage. It is
  set once at the OG seed by the sign of the seed's net vorticity and passed to every
  descendant unchanged (Ch.4, Ch.5): a matter OGU has only matter descendants. So a
  lineage is chirally pure.

  SPIN (the Kerr angular momentum / rotation axis of a given universe) is LOCAL. A
  BHU's spin is the spin of ITS OWN progenitor black hole -- a hole that formed
  astrophysically inside the parent universe, with whatever angular momentum its local
  collapse had. It is NOT the parent universe's net rotation. So spin is re-drawn every
  generation; only the chirality sign is inherited.

Consequences, graded by how true/observable they are:
  * "OGU spin lies in a narrow band" -- TRUE as a derived SHAPE (population.py): viable
    OGU seeds cluster just above the survival threshold |omega/T| > eta_min/C, low-spin
    and exponentially cut off above. Not pinned in absolute units (needs sigma, eta_min).
  * "A BHU's spin follows its parent's spin" -- FALSE: that conflates spin with
    chirality. Only the chirality SIGN is inherited; spin is the local progenitor
    hole's Kerr spin.
  * "spin_us == spin_OGU" -- FALSE: our spin is our progenitor hole's Kerr spin (a
    local astrophysical value, possibly near-extremal), not the OGU's low-band seed spin.
  The genuine observable from OUR spin is a PREFERRED AXIS (our progenitor hole's Kerr
  axis): the CMB "axis of evil" and galaxy-spin handedness should share it (Ch.8).

ETA is NOT up in the air despite the new CMB story. It is anchored TWO independent
ways that agree to ~1%:
  1. CMB: the first/second acoustic peak ratio fixes omega_b -> eta (camb_cmb.py).
  2. BBN: primordial light-element abundances (esp. deuterium D/H) fix eta from
     nuclear physics at z~1e9, long before recombination.
The "CMB = parent Hawking radiation" identification concerns only the ORIGIN of the
thermal plasma at the bounce surface; it does NOT replace recombination, and it leaves
BBN completely untouched. So eta ~ 6.1e-10 stands on two legs, and in SCT it is the
debris dragged through the extrusion from the parent (Ch.4) -- inherited, then read off
identically by both probes.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# --- chirality (inherited) vs spin (local): a one-line truth table ----------
INHERITED = ("chirality sign (matter/antimatter)",)
LOCAL = ("spin magnitude", "spin axis")


def spin_is_inherited() -> bool:
    """Spin is NOT inherited -- only the chirality sign is. (Guards against the
    common conflation.)"""
    return False


def chirality_is_inherited() -> bool:
    """The chirality sign is inherited down a lineage (pure lineages)."""
    return True


# --- eta anchored two ways: CMB (omega_b) and BBN (deuterium) ---------------
ETA10_TO_OMEGABH2 = 1.0 / 273.9     # omega_b h^2 = eta10 / 273.9


def eta10_from_omega_b(omega_b_hsq: float = 0.02237) -> float:
    """Baryon-to-photon ratio in units 1e-10 from the CMB-fixed omega_b."""
    return omega_b_hsq * 273.9


def deuterium_abundance(eta10: float) -> float:
    """Primordial 1e5 * (D/H) from the standard BBN fitting formula
    1e5 D/H = 2.6 (6/eta10)^1.6 (Steigman-type)."""
    return 2.6 * (6.0 / eta10) ** 1.6


def helium_fraction(eta10: float) -> float:
    """Primordial helium mass fraction Y_p (weak eta dependence)."""
    return 0.2470 + 0.0002 * (eta10 - 6.1)


def eta10_from_deuterium(obs_1e5_DH: float = 2.53) -> float:
    """Invert the BBN deuterium relation: eta10 implied by an observed D/H."""
    return brentq(lambda e: deuterium_abundance(e) - obs_1e5_DH, 1.0, 20.0)


def cmb_bbn_concordance(omega_b_hsq: float = 0.02237, obs_1e5_DH: float = 2.53):
    """Return (eta10_CMB, eta10_BBN, fractional difference): the two independent
    determinations of eta and how well they agree."""
    e_cmb = eta10_from_omega_b(omega_b_hsq)
    e_bbn = eta10_from_deuterium(obs_1e5_DH)
    return e_cmb, e_bbn, abs(e_cmb - e_bbn) / e_cmb
