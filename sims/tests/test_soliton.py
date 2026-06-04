"""Tests for the gravity-torsion soliton spectrum and the electroweak S parameter."""

import math

import numpy as np

from cartasis_sims import soliton as so


def test_levels_form_a_rising_tower():
    lv = so.energy_levels(n_levels=5, kind="linear", depth=6.0, E_max=20.0, n_scan=110)
    assert len(lv) >= 4
    assert np.all(np.diff(lv) > 0)           # E_1 < E_2 < ...
    assert np.all(lv > 1.0)                   # above the m=1 rest threshold


def test_linear_confinement_is_regge():
    # a relativistic fermion in a scalar LINEAR well gives E^2 ~ linear in n (Regge)
    lv = so.energy_levels(n_levels=6, kind="linear", depth=6.0, E_max=20.0, n_scan=120)
    slope, intercept, fracrms = so.regge_slope(lv)
    assert fracrms < 0.03                     # E^2 linear to a few percent
    assert slope > 0


def test_anharmonic_well_bends_the_trajectory_up():
    # the bounce (rho^2-wall) well is steeper than linear -> E^2 grows faster than Regge
    lin = so.energy_levels(n_levels=5, kind="linear", depth=6.0, E_max=20.0, n_scan=110)
    bnc = so.energy_levels(n_levels=5, kind="bounce", depth=6.0, E_max=24.0, n_scan=110)
    # compare curvature: bounce E^2 trajectory has larger second difference
    d2_lin = np.diff(lin**2, 2)
    d2_bnc = np.diff(bnc**2, 2)
    assert np.mean(d2_bnc) > np.mean(d2_lin)


def test_effective_mass_profile_is_a_well():
    r = np.linspace(0, 2, 50)
    M = so.effective_mass_profile(r, m=1.0, kind="bounce", depth=6.0)
    assert math.isclose(M[0], 1.0)            # M(0) = m
    assert np.all(np.diff(M) >= -1e-12)       # rises monotonically (a confining wall)


def test_s_parameter_qcd_like_is_in_the_graveyard():
    S = so.s_parameter(so.QCD_FPI_OVER_MV, so.QCD_MA_OVER_MV)
    assert 0.2 < S < 0.32                     # ~0.25, the technicolor graveyard value
    assert S > so.S_LEP_BOUND                  # fails the precision bound


def test_s_below_bound_needs_walking():
    # reaching S < 0.1 requires heavier resonances relative to f_pi (walking)
    need = so.fpi_over_MV_for_S(0.1, so.QCD_MA_OVER_MV)
    assert need < so.QCD_FPI_OVER_MV           # must be smaller than QCD's 0.12
    assert so.s_parameter(need, so.QCD_MA_OVER_MV) < 0.1 + 1e-6
    assert 1.0 / need > 13.0                    # M_V/f_pi > 13 (vs QCD ~8)
