"""Tests for the dark-sector look-back dynamic (gravity felt before light seen)."""

import math

from cartasis_sims import lookback as lb


def test_age_now_is_planck_value():
    assert math.isclose(lb.cosmic_age(1.0) / lb.GYR_S, 13.8, abs_tol=0.3)


def test_particle_horizon_now_is_about_47_gly():
    assert 44.0 < lb.horizon_gly(1.0) < 49.0          # ~46-47 Gly


def test_horizon_grows_with_time():
    # we see further back as the universe ages (the past horizon recedes)
    assert lb.horizon_gly(1.1) > lb.horizon_gly(1.0)


def test_lookback_gain_positive_and_about_1_gly():
    # in +1 Gyr we see ~1 Gly further back -- new parent matter becomes visible
    g = lb.lookback_gain(1.0)
    assert 0.5 < g < 1.5


def test_asymptote_is_finite_and_above_now():
    # de Sitter future: the particle horizon converges to ~64 Gly, not infinity
    h_inf = lb.asymptotic_horizon_gly()
    assert 55.0 < h_inf < 75.0
    assert h_inf > lb.horizon_gly(1.0)
    assert 0.6 < lb.fraction_seen_now() < 0.85         # we see ~74% of it now


def test_event_horizon_is_finite():
    # parent matter beyond the event horizon never becomes visible (only felt)
    eh = lb.event_horizon(1.0) / lb.GLY_M
    assert 10.0 < eh < 25.0                             # ~16-17 Gly comoving
