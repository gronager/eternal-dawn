"""Tests for the primordial spectrum from the bounce (Tier 2 first result)."""

import math

import numpy as np

from cartasis_sims import primordial as pr


def test_p_from_w_matter_and_radiation():
    assert math.isclose(pr.p_from_w(0.0), 2.0)          # matter -> p = 2
    assert math.isclose(pr.p_from_w(1.0 / 3.0), 1.0)    # radiation -> p = 1


def test_analytic_tilt_scale_invariant_at_matter():
    assert math.isclose(pr.analytic_tilt(2.0), 0.0)     # matter: n_T = 0
    assert math.isclose(pr.analytic_tilt(1.0), 2.0)     # radiation: n_T = 2


def test_matter_bounce_is_numerically_scale_invariant():
    # the headline: matter-dominated contraction gives a flat tensor spectrum
    n_T = pr.tilt_numeric(pr.p_from_w(0.0))
    assert abs(n_T) < 0.05


def test_numeric_tilt_matches_analytic():
    for w in (-0.04, 0.0, 0.06, 1.0 / 3.0):
        p = pr.p_from_w(w)
        assert math.isclose(pr.tilt_numeric(p), pr.analytic_tilt(p), abs_tol=0.06)


def test_red_tilt_needs_slightly_soft_contraction():
    # a slightly red spectrum (n_T < 0) requires w < 0 (softer than matter)
    assert pr.tilt_numeric(pr.p_from_w(-0.02)) < 0.0
    assert pr.tilt_numeric(pr.p_from_w(0.02)) > 0.0


def test_power_spectrum_positive():
    ks = np.logspace(-2.0, -0.7, 5)
    P = pr.power_spectrum(ks, pr.p_from_w(0.0))
    assert np.all(P > 0.0)
