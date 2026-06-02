"""Tests for the decisive shared-axis model comparison (SCT vs systematic)."""

import math

import numpy as np

from cartasis_sims import axes as ax


def test_llr_favours_closer_pole():
    # galaxy axis closer to the CMB axis -> positive log-LR (favours SCT)
    assert ax.log_likelihood_ratio(5.0, 30.0, 10.0) > 0
    # closer to the Galactic pole -> negative (favours systematic)
    assert ax.log_likelihood_ratio(30.0, 5.0, 10.0) < 0
    # equidistant -> zero
    assert math.isclose(ax.log_likelihood_ratio(20.0, 20.0, 10.0), 0.0, abs_tol=1e-9)


def test_current_geometry_favours_systematic():
    # galaxy axis ~30 deg from CMB, ~0 from Galactic pole at sigma=5 -> strongly sys
    odds = ax.posterior_odds_sct(30.0, 0.0, 5.0)
    assert odds < 1e-6                          # decisively favours the systematic


def test_sharper_measurement_is_more_decisive():
    # smaller sigma -> odds move further from unity (more decisive)
    o5 = ax.posterior_odds_sct(30.0, 0.0, 5.0)
    o20 = ax.posterior_odds_sct(30.0, 0.0, 20.0)
    assert o5 < o20 < 1.0


def test_sigma_to_decide_for_30deg_poles():
    # poles 30 deg apart need ~10 deg precision for 100:1 odds
    s = ax.sigma_to_decide(30.0, 0.0, 100.0)
    assert 8.0 < s < 12.0
    # closer poles demand sharper data
    assert ax.sigma_to_decide(15.0, 0.0) < ax.sigma_to_decide(45.0, 0.0)


def test_cmb_axis_galactic_pole_separation():
    cmb = ax.lb_to_vec(260.0, 60.0)
    gpole = ax.lb_to_vec(0.0, 90.0)
    sep = ax.discriminating_separation(cmb, gpole)
    assert 25.0 < sep < 35.0                    # ~30 deg apart


def test_symmetric_future_scenario_favours_sct():
    # de-confounded: axis on the CMB axis, 30 deg from Galactic pole, sigma=5
    odds = ax.posterior_odds_sct(0.0, 30.0, 5.0)
    assert odds > 1e6                           # decisively favours SCT
