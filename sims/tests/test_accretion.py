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


def test_mass_budget_branching_caps_and_scales():
    # m = eps f_clean (M_parent/M_vis), capped at the conservation limit
    assert math.isclose(pop.mass_budget_branching(100, 1, 0.1, 1.0), 10.0)
    assert pop.mass_budget_branching(1e9, 1, 0.1, 1.0) <= 1e9   # cap respected
    # equal masses, efficiency 1 -> exactly one viable child (critical)
    assert math.isclose(pop.mass_budget_branching(1, 1, 1.0, 1.0), 1.0)


def test_narrow_band_is_shallow_wide_band_is_deep():
    # narrow band (viable ~ our size, w ~ few) -> m < 1 -> shallow
    m_narrow = pop.mass_budget_branching(3.0, 1.0, 0.1, 1.0)
    assert m_narrow < 1.0
    assert pop.shallow_probability(m_narrow) > 0.8
    # wide band (M_vis << M_parent) -> m >> 1 -> deep
    m_wide = pop.mass_budget_branching(1e4, 1.0, 0.1, 1.0)
    assert m_wide > 1.0
    assert pop.shallow_probability(m_wide) < 1e-3


def test_prob_bhu_stable_for_large_m():
    # log-space normalization: no overflow/nan even for huge branching
    p = pop.prob_bhu(1e5, 1)
    assert 0.0 <= p < 1e-3 and p == p     # finite, not nan


def test_generational_depth_is_geometric_in_epsilon():
    # P(BHU_n|n>=1) = (1-eps) eps^{n-1}: BHU1 dominates, normalized
    n, p = pop.generational_depth_pmf(0.1, n_max=30)
    assert math.isclose(p[0], 0.9, rel_tol=1e-9)         # P(BHU1) = 1 - eps
    assert math.isclose(p[1], 0.09, rel_tol=1e-9)        # P(BHU2)
    assert math.isclose(p.sum(), 1.0, rel_tol=1e-6)


def test_bhu1_independent_of_ogu_size():
    # P(BHU1) depends only on epsilon, NOT on how big the OGU is
    _, p_small = pop.generational_depth_pmf(0.1)
    _, p_big = pop.generational_depth_pmf(0.1)
    assert p_small[0] == p_big[0]                         # OGU size cancelled out


def test_ogu_size_from_sibling_count():
    # M_OGU ~ N_siblings/epsilon * M_vis: 1e11 siblings -> ~1e65 kg for M_vis~1e53
    M_vis = 9.2e52
    M_ogu = 1e11 / 0.1 * M_vis
    assert math.isclose(pop.ogu_siblings(M_ogu, M_vis, epsilon=0.1), 1e11,
                        rel_tol=1e-6)
    assert 1e64 < M_ogu < 1e66                            # ~1e65 kg
