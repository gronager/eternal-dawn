"""Tests for universe-hole growth and the depth it allows."""

import math

from cartasis_sims import growth as g


def test_m_crit_matches_cmb_value_and_scales_inversely():
    # M_crit = hbar c^3/(8 pi G kB T): ~4.5e22 kg at the CMB temperature
    assert math.isclose(g.m_crit(2.7255), 4.5e22, rel_tol=0.1)
    assert g.m_crit(1e12) < g.m_crit(1e9)        # hotter bath -> smaller M_crit


def test_growth_rate_threshold_at_mcrit():
    # dx/dt' = (x^4-1)/x^2: evaporate below x=1, grow above, balance at x=1
    assert g.growth_rate_universal(0.5) < 0
    assert abs(g.growth_rate_universal(1.0)) < 1e-12
    assert g.growth_rate_universal(2.0) > 0


def test_runaway_is_fast_and_finite_above_threshold():
    assert math.isinf(g.runaway_time_universal(0.9))   # below threshold: never
    assert g.runaway_time_universal(2.0) < 1.0          # fast: < tau0
    # for large x0 the blow-up time ~ 1/x0
    assert math.isclose(g.runaway_time_universal(10.0), 0.1, rel_tol=0.05)


def test_depletion_freezes_at_patch_mass():
    # the hole eats its finite local patch and stops at ~ x0 + bath
    t, x, b = g.simulate_depletion(x0=2.0, bath=80.0, t_max=40.0)
    assert math.isclose(x[-1], 82.0, rel_tol=0.05)
    assert b[-1] < 1e-3                                 # bath exhausted


def test_generation_time_scales_as_inverse_T_cubed():
    # tau_gen ~ tau0 ~ t_Hawking(M_crit) ~ M_crit^3 ~ T^-3
    r = g.generation_time(1e12) / g.generation_time(1e13)
    assert math.isclose(r, 1e3, rel_tol=0.05)


def test_causal_rate_integrates_to_horizon_mass():
    # the causal cap c^3/2G times age = the horizon mass: the two models agree
    from cartasis_sims import genesis as gen
    t = 1e6 * gen.HUBBLE_TIME_S
    assert math.isclose(g.causal_growth_rate() * t, gen.horizon_mass(t),
                        rel_tol=1e-9)
    # the rate is ~ a few trillion solar masses per year
    assert 1e12 < g.causal_growth_rate() * 3.156e7 / 1.989e30 < 1e13


def test_max_depth_increases_with_temperature_and_age():
    age = 4.35e17
    assert g.max_depth(age, 1e13) > g.max_depth(age, 1e12)
    assert g.max_depth(10 * age, 1e12) > g.max_depth(age, 1e12)
    # too cold: nothing grows in a Hubble-aged supraverse (sterile foam)
    assert g.max_depth(age, 1e9) == 0
