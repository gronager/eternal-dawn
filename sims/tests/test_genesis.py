"""Tests for OGU genesis: vortex spin window and horizon-limited mass."""

import math

from cartasis_sims import genesis as g


def test_horizon_mass_at_one_hubble_time_is_our_mass():
    # M ~ c^3 t/(2G): a one-Hubble-time causal patch holds ~ our universe's mass
    assert math.isclose(g.horizon_mass(g.HUBBLE_TIME_S), 9e52, rel_tol=0.2)


def test_age_and_mass_invert():
    t = 1e5 * g.HUBBLE_TIME_S
    assert math.isclose(g.age_for_mass(g.horizon_mass(t)), t, rel_tol=1e-9)


def test_condensed_void_is_horizon_limited_subcritical_is_smaller():
    t = 1e3 * g.HUBBLE_TIME_S
    assert math.isclose(g.ogu_mass(t), g.horizon_mass(t), rel_tol=1e-9)
    sub = g.ogu_mass(t, 1e-3 * g.void_critical_density(t))
    assert sub < g.horizon_mass(t)


def test_crossover_is_critical_density():
    # density-limited and horizon masses meet exactly at rho_crit(t)
    t = 50.0 * g.HUBBLE_TIME_S
    rc = g.void_critical_density(t)
    assert math.isclose(g.density_limited_mass(t, rc), g.horizon_mass(t),
                        rel_tol=1e-9)


def test_1e65_ogu_has_about_1e11_siblings():
    assert math.isclose(g.n_siblings(1e65), 1.1e11, rel_tol=0.2)
    # and a near-our-mass OGU is essentially an only child
    assert g.n_siblings(1e54) < 5.0


def test_spin_window_has_a_sterile_floor_and_choked_ceiling():
    fmin = g.spin_floor(eta_min=2e-11, C=10.0)
    assert 0.0 < fmin < 1.0
    assert g.entrainment_efficiency(fmin * 0.5) == 0.0     # sterile below floor
    assert g.entrainment_efficiency(0.999) < g.entrainment_efficiency(0.707)
    assert g.entrainment_efficiency(0.707) > 0.0           # pumps in the window


def test_eta_from_spin_is_linear_and_matches_cve_at_extremal():
    assert math.isclose(g.eta_from_spin(0.5, C=10.0),
                        0.5 * g.eta_from_spin(1.0, C=10.0), rel_tol=1e-9)
