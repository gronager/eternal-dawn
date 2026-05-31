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
