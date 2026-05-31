"""Tests for the void nucleation-and-growth foam (average OGU size)."""

import math

import numpy as np

from cartasis_sims import void_foam as vf


def test_rarer_births_give_bigger_ogus():
    # M_avg ~ (c/beta)^{1/4}: smaller beta -> larger average OGU
    assert vf.avg_ogu_mass(1e-120) > vf.avg_ogu_mass(1e-100)


def test_birth_rate_inverts_average_mass():
    for M in (1e54, 1e60, 1e65):
        b = vf.birth_rate_for_mass(M)
        assert math.isclose(vf.avg_ogu_mass(b), M, rel_tol=1e-6)


def test_fill_scaling_is_quarter_power():
    # t_fill ~ beta^{-1/4}: a 1e4 drop in beta raises t_fill by 10x
    assert math.isclose(vf.fill_time(1e-104) / vf.fill_time(1e-100), 10.0,
                        rel_tol=1e-6)


def test_char_size_is_horizon_of_fill_time():
    b = 1e-110
    # cell radius ~ c t_fill; horizon mass of that radius = c^3 t_fill/2G
    from cartasis_sims import constants as k
    assert math.isclose(vf.char_size(b), k.c * vf.fill_time(b), rel_tol=1e-9)
    assert math.isclose(vf.avg_ogu_mass(b),
                        k.c**2 * vf.char_size(b) / (2 * k.G), rel_tol=1e-9)


def test_johnson_mehl_tessellation_tiles_the_box():
    labels, areas, _ = vf.johnson_mehl_2d(n_seeds=60, grid=120, seed=2)
    assert math.isclose(areas.sum(), 1.0, rel_tol=1e-6)     # cells fill the box
    real = areas[areas > 1e-4]
    assert 5 < len(real) < 60                               # some seeds overgrown
    assert real.std() / real.mean() > 0.2                   # a real size spread
