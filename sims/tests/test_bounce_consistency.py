"""Tests for the bounce-vs-mass-coupling consistency assessment."""

import numpy as np

from cartasis_sims import bounce_consistency as bc


def test_strong_coupling_is_vastly_stronger():
    # the mass-giving coupling (~1/Lambda^2) beats gravitational (~G/c^4) by ~10^33
    r = bc.coupling_ratio()
    assert 1e32 < r < 1e34


def test_strong_coupling_bounce_is_near_nuclear():
    rho = bc.bounce_density("strong")
    # ~0.2x nuclear -- far below the gravitational 1e50
    assert 0.05 < rho / bc.RHO_NUCLEAR < 1.0
    assert rho < bc.bounce_density("gravitational")


def test_neutron_stars_exclude_the_strong_coupling():
    t = bc.strong_coupling_excluded_by_neutron_stars()
    assert t["excluded"]                      # NS core exceeds the strong-bounce density
    assert t["ns_core_exceeds_by"] > 5        # by more than 5x, without bouncing


def test_condensate_melted_at_the_bounce():
    # bounce is far hotter than the condensation scale -> mass sector is radiation there
    assert bc.condensate_is_melted_at_bounce() > 1e3


def test_verdict_names_two_couplings():
    v = bc.summary()["verdict"].lower()
    assert "two couplings" in v and "neutron stars" in v
