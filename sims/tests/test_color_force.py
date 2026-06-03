"""Tests for the colour-channel forces and the confinement gap."""

import math

import numpy as np

from cartasis_sims import color_force as cf


def test_fundamental_casimir():
    assert math.isclose(cf.casimir_fundamental(), 4.0 / 3.0, rel_tol=1e-9)


def test_qcd_color_factors():
    assert math.isclose(cf.color_factor("qqbar_singlet"), -4.0 / 3.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qqbar_octet"), 1.0 / 6.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qq_antitriplet"), -2.0 / 3.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qq_sextet"), 1.0 / 3.0, abs_tol=1e-9)


def test_oge_signs_are_physical():
    # singlet attractive (negative), octet repulsive (positive) at short range
    assert cf.oge_potential(1.0, "qqbar_singlet") < 0
    assert cf.oge_potential(1.0, "qqbar_octet") > 0


def test_only_antisymmetric_channels_bind():
    assert cf.binds_perturbatively("qqbar_singlet")
    assert cf.binds_perturbatively("qq_antitriplet")
    assert not cf.binds_perturbatively("qqbar_octet")
    assert not cf.binds_perturbatively("qq_sextet")


def test_overlap_forces_screen_to_zero():
    # the key honest point: overlap/exchange forces asymptote to ZERO, not ~r
    for ch in cf.CHANNELS:
        assert cf.asymptotes_to_zero(ch)
    r = np.array([1.0, 2.0, 4.0, 8.0])
    Y = np.abs(cf.residual_yukawa(r))
    assert np.all(np.diff(Y) < 0)              # Yukawa screens


def test_confinement_is_not_from_overlap():
    # the linear term RISES (V -> infinity) -- the opposite of the screened overlap;
    # it is imposed for contrast, not produced by the exchange
    r = np.array([1.0, 2.0, 4.0, 8.0])
    assert np.all(np.diff(cf.linear_confinement(r)) > 0)
    # and it dominates the screened pieces at large r
    assert cf.linear_confinement(8.0) > abs(cf.oge_potential(8.0, "qqbar_singlet"))
