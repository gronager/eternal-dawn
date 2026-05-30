"""Tests for the Tier-1 bounce sweep over the equation of state."""

import math

import numpy as np

from cartasis_sims import bounce as bnc


def test_bounce_density_is_universal():
    # a_min = 1 and rho_max = rho_C for every equation of state.
    for w in (0.0, 1.0 / 3.0, 1.0, 5.0 / 3.0):
        sol = bnc.simulate_bounce(w=w, a_init=2.5, t_max=80.0, n_points=12000)
        assert abs(sol.a_min - 1.0) < 1e-3
        assert abs(sol.rho_max - 1.0) < 1e-3


def test_fwhm_matches_exact_scaling():
    # FWHM = 4/n exactly, n = 3(1+w).
    for w in (0.0, 1.0 / 3.0, 1.0):
        sol = bnc.simulate_bounce(w=w, a_init=2.5, t_max=80.0, n_points=12000)
        n = 3.0 * (1.0 + w)
        fwhm = bnc.bounce_fwhm(sol)
        assert math.isclose(fwhm, 4.0 / n, rel_tol=0.01)


def test_max_hubble_is_one_half():
    # max|H| = 1/2 for every w (max of x(1-x) is 1/4).
    for w in (0.0, 1.0 / 3.0, 1.0):
        sol = bnc.simulate_bounce(w=w, a_init=2.5, t_max=80.0, n_points=12000)
        assert math.isclose(np.max(np.abs(sol.H)), 0.5, rel_tol=2e-3)


def test_bounce_sharpens_with_w():
    # stiffer matter -> shorter bounce
    fw = [bnc.bounce_fwhm(bnc.simulate_bounce(w=w, a_init=2.5, t_max=80.0,
                                              n_points=12000))
          for w in (0.0, 1.0 / 3.0, 1.0)]
    assert fw[0] > fw[1] > fw[2]
