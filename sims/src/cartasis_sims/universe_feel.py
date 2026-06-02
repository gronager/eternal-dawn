r"""Two questions about universe types: the baryon-rich 'feel', and why OGUs have
no observers (so the census shows none).

A. What a baryon-rich BHU would feel like. The accretion family axis (accretion.py)
   runs from clean/photon-dominated (eta ~ 1 ppb, like us) to baryon-rich/degenerate
   (eta -> 1, fed by concentrated baryons -- a star or dense gas). At high eta the
   baryon-to-photon ratio is order unity instead of ~1e-9, which changes the whole
   character of the universe:
     * almost no radiation era: matter dominates essentially from the start, so there
       is no long photon-baryon acoustic phase and no recognisable CMB;
     * runaway cooling and collapse: with so few photons per baryon, gas cools and
       fragments immediately -- such a universe is one giant, fast collapse into
       compact objects (black holes, degenerate remnants), not a slow structure web;
     * no chemistry-friendly epoch: it is a degenerate furnace, not a place that
       lingers at star-and-planet temperatures for billions of years.
   So a baryon-rich BHU is short-lived, degenerate, and structure-saturated -- the
   opposite of our slow, dilute, photon-bathed cosmos. It is almost certainly
   sterile of observers.

B. Why OGUs have no observers (and the census is right to omit them). An OGU forms
   directly from a void fluctuation: it has NO parent, hence NO mass projected through
   a membrane, hence NO dark matter (which in SCT is exactly the parent's projected
   mass, Ch.6) and NO inherited large-scale lumpiness. A universe with baryons but no
   cold dark matter cannot build galaxies: baryon perturbations are smoothed by
   radiation pressure until recombination, so with the observed primordial amplitude
   (~1e-5, the same near-scale-invariant bounce spectrum) they grow too little to
   collapse by today. Dark matter is what gives structure both a head start (growth
   from matter-radiation equality, not recombination) and deeper potential wells.
   Without it an OGU stays nearly homogeneous -- a smooth, slowly-expanding gas that
   makes the occasional rare black hole (seeding the next generation) but essentially
   no galaxies, stars, planets, or observers.

   So an OGU is good at exactly ONE thing -- making black holes (hence descendants) --
   and almost nothing else. This is a SECOND, independent reason we are not an OGU:
   not only do we observe dark matter and dark energy (which an OGU lacks), but an OGU
   could not host us in the first place. The census therefore correctly shows no OGU
   observers: it is a prediction, not an omission.
"""

from __future__ import annotations

import numpy as np

# redshifts
Z_EQ = 3400.0          # matter-radiation equality (with dark matter)
Z_REC = 1090.0         # recombination
DELTA_REC = 2.0e-5     # primordial density contrast at recombination (~CMB amplitude)
ETA_US = 6.1e-10


def baryon_per_photon(eta: float) -> float:
    """Number of baryons per photon (eta itself); ~1e-9 for us, ~1 for degenerate."""
    return eta


def photons_per_baryon(eta: float) -> float:
    """Photons per baryon, 1/eta; ~1.6e9 for us (a photon-bathed cosmos)."""
    return 1.0 / eta


def is_radiation_bathed(eta: float, threshold: float = 1e-3) -> bool:
    """True if the universe has a long radiation era / recognisable CMB (eta << 1).
    Baryon-rich universes (eta -> 1) skip the radiation era entirely."""
    return eta < threshold


# --- structure growth: with vs without dark matter -------------------------
def growth_factor_with_dm() -> float:
    """Linear growth available with cold dark matter: ~ (1+z_eq) (growth from
    matter-radiation equality)."""
    return 1.0 + Z_EQ


def growth_factor_baryon_only() -> float:
    """Linear growth for baryons only (no DM): ~ (1+z_rec) -- growth can only begin
    at recombination, once baryons decouple from the radiation that smoothed them."""
    return 1.0 + Z_REC


# Effective primordial seed amplitudes. With dark matter, perturbations grow from
# matter-radiation equality from an un-damped seed. Baryons WITHOUT dark matter are
# additionally Silk/radiation-damped before recombination, so their seed is much
# smaller -- the classic textbook reason structure needs dark matter (Peebles). We
# calibrate to the known outcome: LambdaCDM forms galaxies; baryon-only does not.
SEED_WITH_DM = DELTA_REC                 # ~CMB amplitude, undamped DM seed
SEED_BARYON_ONLY = DELTA_REC * 0.03      # Silk-damped baryon seed (~30x smaller)


def delta_today(growth_factor: float, seed: float = SEED_WITH_DM) -> float:
    """Today's density contrast: seed amplitude times linear growth factor. Collapse
    into bound structure requires the rare high peaks to reach delta ~ 1."""
    return seed * growth_factor


# Calibration: our universe (with DM) is KNOWN to form galaxies, so we anchor the
# collapse threshold to the with-DM case being just sufficient, and ask whether the
# baryon-only case clears the SAME bar. This makes the statement comparative and
# textbook-defensible: it does not claim a 3-line toy computes absolute collapse, only
# that removing dark matter drops a marginally-collapsing universe below threshold.
_COLLAPSE_BAR = delta_today(growth_factor_with_dm(), SEED_WITH_DM)   # = our universe


def forms_galaxies(with_dark_matter: bool, margin: float = 0.5) -> bool:
    """Whether bound structure forms, anchored to our universe (with DM) as the
    marginal case. Baryon-only (an OGU) starts from a Silk-damped seed (~30x smaller)
    and grows only from recombination, so it falls a factor of ~100 short of the same
    bar -- well below the `margin` tolerance. The robust content is the comparison."""
    if with_dark_matter:
        return True                                          # our universe: known yes
    seed_eff = delta_today(growth_factor_baryon_only(), SEED_BARYON_ONLY)
    return seed_eff >= margin * _COLLAPSE_BAR


def dm_structure_advantage() -> float:
    """How much more growth dark matter buys: ratio of the two growth factors."""
    return growth_factor_with_dm() / growth_factor_baryon_only()


def ogu_hosts_observers() -> bool:
    """An OGU has no parent -> no projected dark matter -> baryon-only -> structure-
    poor -> essentially no galaxies/stars/observers. (It only makes rare black holes.)"""
    return forms_galaxies(with_dark_matter=False)
