"""Synthetic-data tests for the lattice analysis extractors (validated independent of any GPU run)."""

import numpy as np

from cartasis_sims import lattice as lat


def test_string_tension_recovered_from_synthetic_potential():
    # build a static potential with a KNOWN string tension and recover it
    r = np.linspace(1.0, 10.0, 20)
    sigma_true, alpha_true, c_true = 0.045, 0.30, 0.55
    V = c_true - alpha_true / r + sigma_true * r
    fit = lat.static_potential_cornell(r, V)
    assert np.isclose(fit["sigma"], sigma_true, rtol=1e-6)
    assert np.isclose(fit["alpha"], alpha_true, rtol=1e-6)
    assert fit["r0_sommer"] > 0          # confining -> a finite Sommer scale


def test_screening_gives_zero_tension():
    # a pure Coulomb (screened) potential has sigma ~ 0 -> no area law
    r = np.linspace(1.0, 10.0, 20)
    V = 0.4 - 0.3 / r
    fit = lat.static_potential_cornell(r, V)
    assert abs(fit["sigma"]) < 1e-6


def test_deconfinement_peak_is_off_grid():
    # a susceptibility peaked between scan points -> beta_c interpolated, not snapped to grid
    beta = np.linspace(5.6, 6.0, 9)
    beta_c_true = 5.785
    chi = 1.0 / (1.0 + ((beta - beta_c_true) / 0.05) ** 2)
    out = lat.deconfinement_beta_c(beta, chi)
    assert np.isclose(out["beta_c"], beta_c_true, atol=0.01)
    assert beta_c_true not in beta        # genuinely between grid points


def test_anomalous_dimension_recovered_from_mode_number():
    # nu(M) = A * M^{4/(1+gamma)} with a known gamma -> recover gamma
    gamma_true = 0.35
    M = np.linspace(0.05, 0.5, 40)
    nu = 12.0 * M ** (4.0 / (1.0 + gamma_true))
    out = lat.anomalous_dimension_from_mode_number(M, nu)
    assert np.isclose(out["gamma_m"], gamma_true, rtol=1e-3)


def test_free_field_mode_number_gives_gamma_zero():
    # a free theory has nu(M) ~ M^4 (gamma = 0)
    M = np.linspace(0.05, 0.5, 40)
    nu = 3.0 * M ** 4
    out = lat.anomalous_dimension_from_mode_number(M, nu)
    assert abs(out["gamma_m"]) < 1e-6


def test_gradient_flow_w0_recovered():
    # construct t^2<E>(t) so that W = t d/dt(t^2 E) crosses 0.3 at a known t -> recover w0
    t = np.linspace(0.01, 4.0, 400)
    # choose t^2 E = 0.15 * t  => W = t * 0.15 = 0.3 at t = 2 => w0 = sqrt(2)
    t2E = 0.15 * t
    w0 = lat.gradient_flow_w0(t, t2E, ref=0.3)
    assert np.isclose(w0, np.sqrt(2.0), rtol=1e-2)
