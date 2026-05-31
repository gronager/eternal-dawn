"""Tests for the accretion-set child-eta map and the BHU-depth read-off."""

import math

from cartasis_sims import accretion as ac
from cartasis_sims import population as pop


def test_fair_sample_is_pure_inheritance():
    # b = 1 (a fair cosmic sample / horizon-scale hole) -> eta_child = eta_parent
    assert math.isclose(ac.eta_child(1.0), ac.ETA_COSMIC, rel_tol=1e-12)


def test_baryon_bias_raises_eta_photon_bias_lowers_it():
    assert ac.eta_child(1e6) > ac.eta_child(1.0) > ac.eta_child(0.1)


def test_eta_is_capped_at_unity():
    # a star (b ~ 1e9) drives a degenerate, baryon-saturated child but eta <= 1
    assert ac.eta_child(1e12) == 1.0
    assert ac.is_baryon_rich(1e9)
    assert not ac.is_baryon_rich(1.0)          # we are in the clean family


def test_our_eta_implies_fair_sample_progenitor():
    # observing eta ~ eta_cosmic means b ~ 1, NOT a star-fed (b~1e9) progenitor
    b_us = ac.ETA_COSMIC / ac.ETA_COSMIC
    assert ac.eta_child(b_us) == ac.ETA_COSMIC
    assert not ac.is_baryon_rich(b_us)


def test_subcritical_branching_is_geometric_and_shallow():
    # m < 1: P(BHU_n | n>=1) = (1-m) m^{n-1}; we are most likely BHU1
    m = 0.4
    assert math.isclose(pop.prob_bhu(m, 1), 1 - m, rel_tol=1e-6)
    assert math.isclose(pop.prob_bhu(m, 2), (1 - m) * m, rel_tol=1e-6)
    assert pop.prob_bhu(m, 1) > pop.prob_bhu(m, 2) > pop.prob_bhu(m, 3)


def test_supercritical_branching_is_deep():
    # m > 1: being BHU1 or BHU2 is exponentially unlikely
    assert pop.shallow_probability(1.8, n_max=2) < 1e-3
    assert pop.shallow_probability(0.2, n_max=2) > 0.9


def test_shallow_probability_decreases_with_m():
    vals = [pop.shallow_probability(m) for m in (0.1, 0.5, 0.9, 1.2)]
    assert vals[0] > vals[1] > vals[2] > vals[3]
