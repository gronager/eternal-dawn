"""Tests for eta anchoring (CMB+BBN) and the spin/chirality distinction."""

import math

from cartasis_sims import spin_eta as se


def test_eta_from_cmb_omega_b():
    # omega_b = 0.02237 -> eta10 ~ 6.1
    assert math.isclose(se.eta10_from_omega_b(0.02237), 6.13, abs_tol=0.1)


def test_cmb_bbn_concordance_to_one_percent():
    e_cmb, e_bbn, diff = se.cmb_bbn_concordance()
    assert 5.8 < e_cmb < 6.4
    assert 5.8 < e_bbn < 6.4
    assert diff < 0.02                         # two independent probes agree to ~1%


def test_deuterium_decreases_with_eta():
    # more baryons -> less leftover deuterium (it gets burned to helium)
    assert se.deuterium_abundance(5.0) > se.deuterium_abundance(7.0)


def test_spin_is_local_chirality_is_inherited():
    # the load-bearing distinction: only the chirality sign is inherited
    assert se.chirality_is_inherited()
    assert not se.spin_is_inherited()


def test_bbn_predicts_observed_deuterium_at_cmb_eta():
    # at the CMB-fixed eta, BBN deuterium matches the observed value within errors
    e_cmb = se.eta10_from_omega_b(0.02237)
    assert abs(se.deuterium_abundance(e_cmb) - 2.53) < 0.1
