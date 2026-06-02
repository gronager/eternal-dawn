"""Tests for the CMB peak-height (baryon-loading) model."""

import math

import numpy as np

from cartasis_sims import cmb_peaks as cp


def test_baryon_loading_is_about_0p6():
    R = cp.baryon_loading(0.02237)
    assert 0.5 < R < 0.75                      # R ~ 0.62 at recombination


def test_first_to_second_peak_ratio_matches_planck():
    # the computed odd/even alternation gives ~2.2 at the observed omega_b
    assert abs(cp.compression_rarefaction_ratio(0.02237) - 2.2) < 0.2


def test_peak_ratio_increases_with_baryon_density():
    assert (cp.compression_rarefaction_ratio(0.028)
            > cp.compression_rarefaction_ratio(0.018))


def test_first_peak_height_near_planck():
    pe, pd = cp.peak_heights()
    # first acoustic peak ~ l 220, D_l ~ 5750 microK^2 (within ~25%)
    assert 180 < pe[0] < 280
    assert 4500 < pd[0] < 7000


def test_finds_five_peaks_in_order():
    pe, pd = cp.peak_heights()
    assert len(pe) == 5
    assert np.all(np.diff(pe) > 0)             # monotonically increasing l


def test_spectrum_is_positive_and_damps():
    ell = np.array([2.0, 220.0, 800.0, 2000.0])
    D = cp.cl_tt(ell)
    assert np.all(D > 0)
    assert D[3] < D[1]                          # Silk-damped tail below the first peak
