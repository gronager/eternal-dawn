"""Tests for the Fierz projection of the torsion four-fermion interaction."""

import math

from cartasis_sims import fierz as fz


def test_code_validations():
    # the V-A self-Fierz identity and Fierz^2 = identity must hold exactly
    assert fz.validate_VmA_self_fierz() < 1e-9
    assert fz.validate_involution() < 1e-9


def test_hehl_datta_channel_couplings():
    c = fz.hehl_datta_channels()
    assert math.isclose(c["S"], 0.25, abs_tol=1e-9)
    assert math.isclose(c["P"], -0.25, abs_tol=1e-9)
    assert math.isclose(c["V"], 0.5, abs_tol=1e-9)
    assert math.isclose(c["A"], 0.5, abs_tol=1e-9)
    assert abs(c["T"]) < 1e-9


def test_vector_axial_couplings_are_equal():
    # the decisive, walking-favourable result: G_A/G_V = 1 (vs QCD-like != 1)
    assert math.isclose(fz.GA_over_GV(), 1.0, abs_tol=1e-9)


def test_scalar_channel_is_attractive():
    # G_S > 0 drives the chiral condensate the soliton forms
    assert fz.hehl_datta_channels()["S"] > 0
