"""Tests for the acoustic-peak scaffold (sound horizon, peak positions)."""

import numpy as np

from cartasis_sims import acoustic as ac


def test_sound_horizon_is_about_145_mpc():
    rs = ac.sound_horizon()
    assert 135.0 < rs < 155.0                 # Planck ~ 144.4 Mpc


def test_distance_to_last_scattering():
    dc = ac.comoving_distance_lss()
    assert 13000.0 < dc < 14500.0             # Planck ~ 13870 Mpc


def test_acoustic_scale_matches_planck():
    theta, lA, rs, dc = ac.acoustic_scale()
    # theta_s ~ 0.0104 rad; acoustic multipole l_A ~ 300
    assert abs(theta - 0.0104) < 0.0006
    assert 280.0 < lA < 320.0


def test_peak_positions_match_planck_within_few_percent():
    sct = ac.peak_multipoles(n_peaks=5)
    planck = np.array([220.6, 537.5, 810.8, 1120.9, 1444.0])
    # first peak near 220; all within ~6%
    assert abs(sct[0] - 220.0) < 25.0
    assert np.all(np.abs(sct - planck) / planck < 0.07)


def test_schematic_cl_peaks_sit_at_computed_multipoles():
    ell = np.arange(2, 2000.0)
    Dl = ac.schematic_cl(ell)
    assert np.all(Dl >= 0.0)
    # the dominant (first) peak is within a bin of the first acoustic multipole
    l1 = ac.peak_multipoles(n_peaks=1)[0]
    peak_ell = ell[200:1000][np.argmax(Dl[200:1000])]
    assert abs(peak_ell - l1) < 40.0
