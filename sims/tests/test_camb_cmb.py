"""Tests for the CAMB precision CMB spectrum (skipped if CAMB is not installed)."""

import math

import numpy as np
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


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_dark_matter_sets_peak_heights():
    # more cold dark matter -> less radiation driving -> lower first peak
    p_lo = cc.peak_table(omch2=0.09)
    p_hi = cc.peak_table(omch2=0.15)
    assert p_lo[0][1] > p_hi[0][1]                 # first peak falls with omega_c
    # and the first/third ratio drops (3rd lifts relative to 1st)
    assert p_lo[0][1] / p_lo[2][1] > p_hi[0][1] / p_hi[2][1]


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_dark_energy_sets_peak_positions():
    # more dark energy (higher H0 at fixed physical densities) slides peaks to lower l
    l1_lo = cc.peak_table(H0=60.0)[0][0]
    l1_hi = cc.peak_table(H0=75.0)[0][0]
    assert l1_hi < l1_lo                            # peak position shifts with Omega_L


def test_planck_binned_tt_loads():
    # the real Planck 2018 binned TT band powers ship in data/planck/
    ell, Dl, sigma, best = cc.planck_binned_tt()
    assert ell.size > 70                            # ~83 band powers
    assert ell[0] < 60 and ell[-1] > 2000          # l ~ 48 .. 2500
    assert np.all(sigma > 0) and np.all(Dl > 0)
    # first acoustic peak shows up in the data near l~220, D~5800 microK^2
    j = int(np.argmin(np.abs(ell - 220)))
    assert 5000 < Dl[j] < 6500


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_ed_prediction_fits_real_planck():
    # ED-sourced parameters (not re-fit) reproduce the real binned spectrum at chi2/N~1
    chi2, n = cc.chi2_vs_planck()
    assert n > 70
    assert chi2 / n < 1.3                           # consistency hurdle cleared


@pytest.mark.skipif(not cc.HAVE_CAMB, reason="CAMB not installed")
def test_peak_height_ratios_match_planck():
    r = cc.peak_height_ratios()
    assert abs(r["r12"] - 2.2) < 0.3               # odd/even (baryon) ratio ~2.2
    assert abs(r["r13"] - 2.3) < 0.4               # 1st/3rd (dark-matter) ratio ~2.3


def test_baryon_fraction_is_about_one_sixth():
    # omega_c / omega_b = 5.36 -> baryon fraction f ~ 0.157 ~ 1/6 (the target Ch.6 owes)
    ob, oc = cc.SCT_PARAMS["ombh2"], cc.SCT_PARAMS["omch2"]
    f = ob / (ob + oc)
    assert abs(f - 1.0 / 6.0) < 0.02
