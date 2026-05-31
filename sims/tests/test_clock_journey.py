"""Tests for the clock-journey timescales."""

from cartasis_sims import os_collapse as osc
from cartasis_sims import constants as k


def test_infall_proper_time_scales():
    # stellar 10 Msun ~ 0.15 ms; universe-mass ~ tens of Gyr
    t_stellar = osc.infall_proper_time(10 * k.M_sun)
    assert 1e-4 < t_stellar < 1e-3
    t_uni = osc.infall_proper_time(9.246e52) / k.year
    assert 1e10 < t_uni < 1e11          # ~22 Gyr


def test_infall_linear_in_mass():
    assert abs(osc.infall_proper_time(2e31) / osc.infall_proper_time(1e31) - 2.0) < 1e-9


def test_hawking_cubic_in_mass_and_huge():
    th1 = osc.hawking_time(1e31)
    th2 = osc.hawking_time(2e31)
    assert abs(th2 / th1 - 8.0) < 1e-6   # M^3
    # 10 Msun evaporation ~ 1e70 yr
    assert 1e69 < osc.hawking_time(10 * k.M_sun) / k.year < 1e71


def test_return_lag_factor_scales_as_M_squared():
    # Hawking/infall ratio ~ M^2
    r1 = osc.hawking_time(1e31) / osc.infall_proper_time(1e31)
    r2 = osc.hawking_time(1e32) / osc.infall_proper_time(1e32)
    assert abs(r2 / r1 - 100.0) < 1e-3
