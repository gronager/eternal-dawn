"""Tests for the void weight and the de Sitter cap on OGU growth."""

import math

from cartasis_sims import void_scale as vs


def test_vacuum_density_is_dark_energy_scale():
    # rho_Lambda ~ 6e-27 kg/m^3 (a few proton masses per cubic metre): nonzero weight
    rho = vs.vacuum_density()
    assert 4e-27 < rho < 8e-27
    assert 1.0 < vs.vacuum_energy_gev_m3() < 10.0     # ~3 GeV/m^3


def test_void_has_nonzero_weight_so_not_polygons():
    # the whole point: nonzero rho_Lambda -> nonzero gravity-scaled volume -> gaps
    assert vs.vacuum_density() > 0.0


def test_nariai_mass_caps_ogu_near_hubble_mass():
    # the de Sitter Nariai bound caps OGU growth at ~ our mass (~1e53 kg)
    M = vs.nariai_mass()
    assert 1e52 < M < 1e53
    # Nariai is below the de Sitter horizon mass
    assert M < vs.de_sitter_horizon_mass()


def test_void_dominates_in_the_dilute_regime():
    # far-apart islands -> void holds essentially all the gravitational content
    assert vs.gravity_scaled_void_fraction(1e6) > 0.999999
    # touching (separation ~1) -> void and universes comparable
    assert 0.3 < vs.gravity_scaled_void_fraction(1.0) < 0.7
