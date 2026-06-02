r"""Child eta as a function of what the parent black hole eats.

The extruder result (extruder.py) is that the bounce is adiabatic, D = 1: entropy
is conserved across the crossing. So the child does NOT regenerate its photons at
the bounce temperature -- it inherits the baryon-to-photon (baryon-to-entropy)
ratio of the MATERIAL THE HOLE ACTUALLY ATE:

    eta_child = eta_infall = b * eta_parent,

where b = eta_infall / eta_cosmic is the baryon enrichment of the accreted
material relative to a fair cosmic sample of the parent:

  * b = 1  -- a fair sample (incl. the photon bath): pure inheritance,
              eta_child = eta_parent. A horizon/Hubble-scale hole contains a fair
              sample by construction; a radiation-era primordial hole captures
              plasma whose eta is the same comoving eta. Both give b ~ 1.
  * b > 1  -- baryon-biased accretion (stars, gas, galaxies have few thermal
              photons per baryon): a MORE matter-rich child, eta_child up to ~1
              (a fully baryon-dominated, degenerate universe).
  * b < 1  -- photon-biased accretion (void / intergalactic medium, radiation):
              a CLEANER, more photon-dominated child, eta_child < eta_parent.

eta is capped at ~1 (one net baryon per quantum: a degenerate, baryon-saturated
child). This is exactly the user's point that a growing hole "pulls in stuff post
the primordial annihilation, so there is another eta inside than how the parent
was born" -- made quantitative, and signed: which way eta drifts depends on the
accretion bias.

Two families follow. A star/galaxy-fed (baryon-biased) hole spawns a baryon-rich
/ degenerate universe (eta ~ 1); a fair-sample or radiation/void-fed hole spawns a
clean, photon-dominated universe (eta << 1). OUR universe is photon-dominated
(eta ~ 6e-10 ~ eta_cosmic), so our progenitor hole was NOT a baryon-biased stellar
accretor: it grew by fair-sample or radiation-dominated capture -- a
horizon-scale or primordial hole. That is a structural prediction about our
parentage, not a free choice.

Contrast (non-adiabatic): had the bounce instead fully re-thermalized the infall
to T_bounce (D >> 1), the child would carry eta ~ b_rest * T_bounce / (m_p c^2)
with b_rest the rest-mass baryon loading -- a different, regeneration-set value.
The adiabatic result is what makes eta_child = eta_infall; see extruder.py.
"""

from __future__ import annotations

import math

from . import cve_filter as cve

ETA_COSMIC = cve.ETA_OBS              # our parent's (and our) cosmic baryon-to-photon ratio
T_BOUNCE_GEV = 1.053e7               # plasma temperature at rho_C ~ 1e50 (cve_filter)
M_PROTON_GEV = 0.938

# Representative accretion sources, as enrichment b = eta_infall / eta_cosmic.
# Order-of-magnitude; the point is the ordering and which side of b=1 we sit on.
SOURCES = {
    "void / IGM (photon-biased)": 0.1,
    "radiation-era plasma":       1.0,
    "fair cosmic sample":         1.0,
    "horizon-scale hole":         1.0,
    "galaxy ISM (baryon-biased)": 3.0e6,
    "gas cloud":                  1.0e8,
    "a star":                     1.0e9,
}


def eta_child(b: float, eta_parent: float = ETA_COSMIC) -> float:
    """Adiabatic inheritance: eta_child = min(b * eta_parent, 1)."""
    return min(b * eta_parent, 1.0)


def is_baryon_rich(b: float, eta_parent: float = ETA_COSMIC,
                   threshold: float = 1.0e-2) -> bool:
    """Child is in the baryon-rich/degenerate family if eta_child is driven up
    near order unity (a baryon-saturated, degenerate universe). Only concentrated
    baryon feeding (b ~ 1e8-1e9, a gas cloud or star) reaches this."""
    return eta_child(b, eta_parent) > threshold


def eta_child_thermalized(b_rest: float,
                          T_bounce_GeV: float = T_BOUNCE_GEV) -> float:
    """Non-adiabatic contrast: if the bounce re-thermalized the infall to
    T_bounce, eta ~ b_rest * T_bounce / (m_p c^2). Shown only to make clear that
    the adiabatic D=1 result is what gives eta_child = eta_infall instead."""
    return min(b_rest * T_bounce_GeV / M_PROTON_GEV, 1.0)
