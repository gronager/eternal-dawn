"""Tests for the CVE asymmetry estimate."""

import math

from cartasis_sims import cve_filter as cve


def test_bounce_temperature_is_about_1e7_GeV():
    T = cve.bounce_temperature_GeV(1.0e50)
    assert 1e6 < T < 1e8          # ~10^7 GeV: EW restored, sphalerons active


def test_vorticity_over_T_is_about_1e_minus_11():
    x = cve.vorticity_over_T(1.0e50, spin_fraction=1.0)
    assert 1e-12 < x < 1e-10


def test_asymmetry_scales_with_C_and_spin():
    base = cve.asymmetry_estimate(spin_fraction=0.5, C=1.0)
    assert math.isclose(cve.asymmetry_estimate(spin_fraction=0.5, C=2.0),
                        2.0 * base, rel_tol=1e-9)
    assert math.isclose(cve.asymmetry_estimate(spin_fraction=1.0, C=1.0),
                        2.0 * base, rel_tol=1e-9)


def test_needed_factor_is_order_tens():
    # eta_obs is ~1-2 orders above omega/T, so the needed C*spin is ~tens
    need = cve.spinC_needed(1.0e50)
    assert 5 < need < 500


def test_adiabatic_inheritance_reproduces_parent():
    # D=1 (entropy-conserving extrusion) passes the parent's ratio straight through
    assert math.isclose(cve.inherited_eta(cve.ETA_OBS, dilution=1.0),
                        cve.ETA_OBS, rel_tol=1e-12)


def test_dilution_suppresses_inherited_eta():
    # larger entropy multiplication drives eta down the lineage ("hell downward")
    assert cve.inherited_eta(cve.ETA_OBS, dilution=1e3) < cve.ETA_OBS
    assert math.isclose(cve.inherited_eta(cve.ETA_OBS, dilution=1e3),
                        cve.ETA_OBS / 1e3, rel_tol=1e-12)


def test_horizon_dilution_is_catastrophic():
    # if the interior thermalized to the parent horizon entropy, D is enormous
    # (~1e49 for a 1e53 kg parent) and inherited eta is crushed far below eta_obs;
    # since we DO observe eta~6e-10, the bounce must be near-adiabatic.
    D = cve.dilution_horizon(1.0e53)
    assert D > 1e30
    assert cve.inherited_eta(cve.ETA_OBS, dilution=D) < 1e-40


def test_combined_eta_takes_the_winner():
    # inheritance dominates at low D; fresh skew floors it at high D
    fresh = cve.asymmetry_estimate(spin_fraction=1.0, C=1.0)
    assert math.isclose(cve.combined_eta(dilution=1.0), cve.ETA_OBS, rel_tol=1e-9)
    assert math.isclose(cve.combined_eta(dilution=1e40), fresh, rel_tol=1e-9)
