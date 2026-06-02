"""Tests for the void max-entropy eigenstate and OGUs as suppressed fluctuations."""

import math

from cartasis_sims import void_eigenstate as ve


def test_de_sitter_entropy_is_about_1e122():
    # the Gibbons-Hawking entropy of our de Sitter horizon ~ 1e122 k_B
    S = ve.de_sitter_horizon_entropy()
    assert 1e121 < S < 1e123


def test_empty_de_sitter_is_the_maximum():
    # inserting any sub-Nariai mass LOWERS the total horizon entropy
    assert ve.empty_is_maximum()
    assert ve.entropy_deficit(0.01) > 0
    assert ve.entropy_deficit(0.5) > 0


def test_deficit_grows_with_mass():
    assert ve.entropy_deficit(0.5) > ve.entropy_deficit(0.05)


def test_small_mass_deficit_matches_boltzmann_exponent():
    # the entropy deficit -> 2 pi mu (the de Sitter Boltzmann exponent = instanton I)
    assert math.isclose(ve.entropy_deficit(1e-3), ve.deficit_small_mass(1e-3),
                        rel_tol=1e-2)
    assert ve.nucleation_weight_is_entropy_deficit()


def test_deficit_is_a_small_fraction_of_total_entropy():
    # even a Nariai-mass insertion is a modest fraction of S_dS (the eigenstate is
    # only slightly perturbed) -- the deficit at small mass is tiny
    assert ve.entropy_deficit(1e-4) < 1e-3
