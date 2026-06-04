"""Tests for the leading-order electroweak S parameter from the soliton sector."""

import math

import numpy as np

from cartasis_sims import electroweak_S as ew


def test_spectral_difference_is_positive_and_finite():
    M = 1.0
    s = np.linspace(4 * M**2 * 1.001, 300 * M**2, 3000)
    d = ew.spectral_difference(s, M)
    assert np.all(d >= 0)                         # rho_V - rho_A > 0 -> S > 0
    assert np.isfinite(np.trapezoid(d / s**2, s))  # UV-finite (chiral-odd)


def test_vector_exceeds_axial_near_threshold():
    M = 1.0
    s = np.array([5.0, 10.0, 30.0]) * M**2
    assert np.all(ew.spectral_V(s, M) > ew.spectral_A(s, M))


def test_spectral_difference_vanishes_at_high_s():
    # chiral symmetry restored at high s: the difference ~ M^2/s -> 0
    M = 1.0
    assert ew.spectral_difference(np.array([4.5]) * M**2, M)[0] > \
        ew.spectral_difference(np.array([400.0]) * M**2, M)[0]


def test_leading_S_is_in_the_graveyard():
    assert math.isclose(ew.s_leading(3, 1), 3 / (6 * math.pi), rel_tol=1e-9)
    v = ew.verdict(3, 1)
    assert v["in_graveyard"]                      # > 0.1
    assert v["S_leading"] > ew.S_LEP_BOUND


def test_more_colour_makes_S_worse():
    # Pati-Salam (N_c=4) raises S -- the extra colour is unfavourable
    assert ew.s_leading(4, 1) > ew.s_leading(3, 1)


def test_walking_can_reach_the_bound():
    # S < 0.1 requires M_V/f_pi >~ 14 (a ~1.7x walk over QCD's ~8)
    need = ew.verdict(3, 1)["MV_over_fpi_needed"]
    assert need > 13.0
    assert ew.s_walking(1.0 / need) < 0.1 + 1e-6
