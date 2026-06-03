"""Tests for the cosmic-birefringence EB/TB fingerprint."""

import math

import numpy as np
import pytest

from cartasis_sims import birefringence as bi
from cartasis_sims import camb_cmb as cc


def test_no_rotation_no_parity_odd_power():
    # LambdaCDM (beta = 0): EB and TB vanish exactly
    ee = np.array([10.0, 40.0, 5.0])
    bb = np.array([0.1, 0.5, 0.2])
    te = np.array([100.0, -50.0, 30.0])
    assert np.allclose(bi.eb_spectrum(ee, bb, 0.0), 0.0)
    assert np.allclose(bi.tb_spectrum(te, 0.0), 0.0)


def test_eb_matches_half_sin4beta_ee():
    ee = np.array([42.0])
    bb = np.array([0.0])
    beta = 0.3
    expected = 0.5 * math.sin(4 * math.radians(beta)) * 42.0
    assert math.isclose(bi.eb_spectrum(ee, bb, beta)[0], expected, rel_tol=1e-12)


def test_tb_matches_sin2beta_te():
    te = np.array([-130.0])
    beta = 0.3
    expected = math.sin(2 * math.radians(beta)) * (-130.0)
    assert math.isclose(bi.tb_spectrum(te, beta)[0], expected, rel_tol=1e-12)


def test_amplitude_grows_with_beta():
    ee = np.array([42.0]); bb = np.array([0.0])
    assert (bi.eb_peak_amplitude(ee, bb, 0.5)
            > bi.eb_peak_amplitude(ee, bb, 0.3)
            > bi.eb_peak_amplitude(ee, bb, 0.1))


def test_combined_beta_is_about_a_third_of_a_degree():
    beta, sigma, _ = bi.combined_beta()
    assert 0.25 < beta < 0.40            # ~0.3 deg consensus
    assert sigma > 0


def test_axis_map_dipole_aligned_with_axis():
    lon, lat, beta = bi.axis_anisotropy_map(0.3, 0.2, axis_lonlat=(260.0, 60.0))
    # the maximum of the dipole field sits toward the axis direction (b=60)
    j, i = np.unravel_index(np.argmax(beta), beta.shape)
    assert lat[j] > 30.0                 # peak in the northern axis hemisphere


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_real_spectra_give_expected_eb_amplitude():
    ell, TT, EE, BB, TE = cc.cmb_polarization_spectra(lmax=2000)
    # beta ~ 0.3 deg leaves a sub-microK but nonzero EB peak from the real EE
    amp = bi.eb_peak_amplitude(EE, BB, 0.3)
    assert 0.2 < amp < 1.0
