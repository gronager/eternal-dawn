"""Tests for the timescape clock-rate resolution of the Hubble tension."""

import math

from cartasis_sims import hubble_timescape as ht


def test_homogeneous_limit_has_no_effect():
    # R(0) = 1: with no voids the dressed and bare Hubble coincide
    assert math.isclose(ht.dressed_over_bare(0.0), 1.0, rel_tol=1e-12)


def test_dressed_over_bare_at_bestfit_void_fraction():
    # Wiltshire closed form at f_v0 = 0.762 -> ~1.28 (a 28% wall-vs-average excess)
    assert math.isclose(ht.dressed_over_bare(0.762), 1.2825, abs_tol=1e-3)


def test_implied_bare_hubble_matches_literature():
    # dressed 61.7 / R(0.762) -> bare ~48, the published timescape bare H0
    bare = ht.H0_DRESSED_BESTFIT / ht.dressed_over_bare(0.762)
    assert 47.0 < bare < 50.0


def test_effect_has_right_sign_and_enough_magnitude():
    # local/wall exceeds global/average, by more than the observed gap
    assert ht.wall_vs_average_offset() > ht.observed_tension_fraction() > 0.0
    assert ht.resolves_in_right_direction()


def test_modest_void_fraction_already_matches_the_tension():
    # the void fraction that reproduces the exact 8.4% gap is ~0.47, BELOW the
    # observed 0.76 -- the real void fraction has headroom
    f_needed = ht.void_fraction_to_match_tension()
    assert 0.4 < f_needed < 0.55
    assert f_needed < ht.FV0_BESTFIT


def test_early_universe_is_untouched():
    # f_v -> 0 at recombination, so the CMB sound-horizon H0 is unmodified
    assert ht.void_fraction_at_recombination() < 1e-3
    assert ht.early_universe_unaffected()


def test_void_fraction_tracker_normalisation():
    # f_v(t0) = f_v0 and f_v -> 0 as t -> 0
    assert math.isclose(ht.void_fraction(1.0), ht.FV0_BESTFIT, rel_tol=1e-9)
    assert ht.void_fraction(1e-6) < 1e-5
