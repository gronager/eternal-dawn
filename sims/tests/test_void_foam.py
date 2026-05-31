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


def test_ogus_are_effectively_immortal():
    # Hawking lifetime of OGU masses dwarfs any supraverse timescale -> no
    # evaporation death; OGUs only coarsen by merging.
    assert vf.ogu_evaporation_lifetime(9.2e52) > 1e140      # our-mass: ~1e143 s
    assert vf.ogu_evaporation_lifetime(1e65) > 1e170        # ~1e178 s
    # scales as M^3
    assert math.isclose(vf.ogu_evaporation_lifetime(2e53)
                        / vf.ogu_evaporation_lifetime(1e53), 8.0, rel_tol=1e-6)


def test_foam_bath_mass_equal_masses():
    # equal masses -> bath equivalent mass equals them (all at marginal balance)
    assert math.isclose(vf.foam_bath_mass([1.0, 1.0, 1.0, 1.0]), 1.0, rel_tol=1e-9)


def test_coarsening_grows_big_evaporates_small():
    # negative-heat-capacity instability: above M_bar grow, below shrink
    masses = np.array([0.7, 1.0, 1.5])
    Mb = vf.foam_bath_mass(masses)
    rate = vf.coarsening_rate(masses)
    assert rate[masses > Mb][0] > 0.0          # big OGU grows
    assert rate[masses < Mb][0] < 0.0          # small OGU evaporates


def test_coarsening_conserves_mass_and_coarsens():
    rng = np.random.default_rng(0)
    m0 = rng.lognormal(0.0, 0.25, 12)
    m0 /= m0.mean()
    tau, M = vf.simulate_coarsening(m0, tau_max=4.0)
    alive0, alive1 = (M[:, 0] > 0.02).sum(), (M[:, -1] > 0.05).sum()
    assert alive1 < alive0                      # number drops (some die)
    assert M[:, -1].max() > M[:, 0].max()       # the biggest grows
    # surviving foam conserves mass to within the small dead residuals
    assert abs(M[:, -1].sum() - M[:, 0].sum()) < 0.3


def test_johnson_mehl_tessellation_tiles_the_box():
    labels, areas, _ = vf.johnson_mehl_2d(n_seeds=60, grid=120, seed=2)
    assert math.isclose(areas.sum(), 1.0, rel_tol=1e-6)     # cells fill the box
    real = areas[areas > 1e-4]
    assert 5 < len(real) < 60                               # some seeds overgrown
    assert real.std() / real.mean() > 0.2                   # a real size spread
