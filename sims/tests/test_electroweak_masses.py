"""Tests for the W, Z, Higgs masses in the composite/walking reading."""

import math

from cartasis_sims import electroweak_masses as ew


def test_W_mass_within_one_percent():
    s = ew.summary()["m_W"]
    assert s["off"] < 1.01            # m_W = (1/2) g v lands within ~1%


def test_Z_mass_within_one_percent():
    s = ew.summary()["m_Z"]
    assert s["off"] < 1.01


def test_custodial_rho_is_one():
    # the clean structural prediction of a doublet condensate
    assert math.isclose(ew.rho_parameter(), 1.0, rel_tol=1e-6)


def test_custodial_relation_recovers_mW():
    # m_W = m_Z cos theta_W
    assert math.isclose(ew.m_W_from_custodial(), ew.m_W(), rel_tol=1e-6)


def test_observed_higgs_sits_at_the_walking_edge():
    lo, hi = ew.higgs_bracket()
    # observed Higgs is inside the bracket and near the light (walking-dilaton) end
    assert lo * 0.95 < ew.OBSERVED["H"] < hi
    assert ew.OBSERVED["H"] < 0.5 * hi          # far below the heavy-sigma end


def test_higgs_ratio_is_about_half_v():
    assert math.isclose(ew.OBSERVED["H"] / ew.V_EW, 0.51, abs_tol=0.02)
