"""Tests for the baryon-rich feel and OGU structure-poverty."""

import math

from cartasis_sims import universe_feel as uf


def test_we_are_photon_bathed_degenerate_is_not():
    assert uf.is_radiation_bathed(uf.ETA_US)             # eta ~ 1ppb: long radiation era
    assert not uf.is_radiation_bathed(0.6)               # eta -> 1: no radiation era
    # photons per baryon: ~1.6e9 for us, ~order 1 for degenerate
    assert uf.photons_per_baryon(uf.ETA_US) > 1e8
    assert uf.photons_per_baryon(0.6) < 10


def test_dark_matter_gives_a_growth_advantage():
    assert uf.growth_factor_with_dm() > uf.growth_factor_baryon_only()
    assert uf.dm_structure_advantage() > 2.0            # ~3x head start


def test_bhu_forms_galaxies_ogu_does_not():
    assert uf.forms_galaxies(with_dark_matter=True)      # our universe (anchored yes)
    assert not uf.forms_galaxies(with_dark_matter=False)  # OGU: baryon-only, falls short


def test_ogu_hosts_no_observers():
    # the headline: an OGU has no dark matter -> no galaxies -> no observers
    assert not uf.ogu_hosts_observers()


def test_ogu_far_below_collapse_bar():
    # OGU clears only a small fraction of our universe's collapse threshold
    frac = (uf.delta_today(uf.growth_factor_baryon_only(), uf.SEED_BARYON_ONLY)
            / uf._COLLAPSE_BAR)
    assert frac < 0.1                                    # ~1% -- nowhere near collapse
