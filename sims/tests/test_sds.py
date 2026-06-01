"""Tests for Schwarzschild-de Sitter two-horizon thermodynamics."""

import math

import numpy as np

from cartasis_sims import sds


def test_two_horizons_below_nariai():
    h = sds.sds_horizons(0.5 * sds.MU_NARIAI)
    assert h is not None
    rb, rc = h
    assert 0 < rb < rc                          # BH horizon inside cosmological one


def test_horizons_merge_at_nariai():
    h = sds.sds_horizons(sds.MU_NARIAI * (1 - 1e-9))
    rb, rc = h
    assert math.isclose(rb, sds.R_NARIAI, rel_tol=1e-3)
    assert math.isclose(rc, sds.R_NARIAI, rel_tol=1e-3)


def test_no_horizons_above_nariai():
    assert sds.sds_horizons(1.2 * sds.MU_NARIAI) is None


def test_bh_always_hotter_than_cosmological_horizon():
    # T_b > T_c for every sub-Nariai mass -> net outflow -> evaporation
    for f in (0.05, 0.3, 0.6, 0.9, 0.999):
        mu = f * sds.MU_NARIAI
        assert sds.temperature_ratio(mu) > 1.0
        assert sds.net_outflow(mu)


def test_true_equilibrium_only_at_nariai():
    assert not sds.in_true_equilibrium(0.9 * sds.MU_NARIAI)
    # only the exact Nariai coincidence gives T_b = T_c
    assert sds.in_true_equilibrium(sds.MU_NARIAI, tol=1e-3)


def test_temperature_gap_closes_toward_nariai():
    # near Nariai T_b/T_c -> 1 (lukewarm plateau that mimics equilibrium)
    assert (sds.temperature_ratio(0.999 * sds.MU_NARIAI)
            < sds.temperature_ratio(0.5 * sds.MU_NARIAI))
    assert sds.temperature_ratio(0.99999 * sds.MU_NARIAI) < 1.01
