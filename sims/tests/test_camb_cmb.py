"""Tests for the CAMB precision CMB spectrum (skipped if CAMB is not installed)."""

import math

import pytest

from cartasis_sims import camb_cmb as cc

# the eta -> omega_b cross-check needs no CAMB; the spectrum tests do.
pytestmark_camb = pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")


def test_eta_predicts_planck_baryon_density():
    # inherited eta ~ 6.1e-10 -> omega_b ~ 0.0223, within 0.5% of Planck's 0.02237
    omb = cc.baryon_density_from_eta(6.1e-10)
    assert math.isclose(omb, 0.02237, rel_tol=0.01)


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_first_peak_position_and_height():
    pe = cc.peak_table()
    l1, d1 = pe[0]
    assert abs(l1 - 220) < 15                      # first acoustic peak ~ 220
    assert 5400 < d1 < 6000                         # ~5700 microK^2


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_five_peaks_track_planck():
    pe = cc.peak_table()
    planck = [(220, 5710), (540, 2560), (810, 2480), (1130, 1240), (1450, 800)]
    assert len(pe) == 5
    for (lm, dm), (lp, dp) in zip(pe, planck):
        assert abs(lm - lp) < 35                    # positions within a few %
        assert abs(dm - dp) / dp < 0.10             # heights within ~10%


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_spectrum_is_finite_and_damped():
    import numpy as np
    ell, TT = cc.cmb_tt_spectrum(lmax=2200)
    assert np.all(np.isfinite(TT))
    # Silk-damped tail: power at l~2000 well below the first peak
    assert TT[2000] < 0.2 * TT[220]
