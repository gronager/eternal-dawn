"""Tests for the scale-gap dimensional-transmutation resolution."""

import numpy as np

from cartasis_sims import scale_gap as sg


def test_transmutation_inverts_consistently():
    # Lambda(g2) and g2_for_scale(Lambda) are inverses
    g2 = sg.g2_for_scale(Lambda=173.0, b0=7.0)
    assert np.isclose(sg.transmuted_scale(g2, b0=7.0), 173.0, rtol=1e-6)


def test_required_coupling_is_ordinary():
    # the 10^17 hierarchy needs only an O(0.01-0.2) coupling, not a tuning
    for b0, d in sg.hierarchy_is_ordinary().items():
        assert 0.005 < d["g_T^2"] < 0.3
        assert d["alpha_T"] < 0.05


def test_coupling_is_asymptotically_free():
    # weak in the UV (M_Pl), strong in the IR (Lambda)
    a_lambda = sg.running_coupling(173.0)
    a_planck = sg.running_coupling(sg.M_PL_GEV)
    assert a_lambda > a_planck
    assert a_planck < 0.1            # few-percent at the Planck scale


def test_rescaled_top_lands_in_ballpark():
    s = sg.scale_summary()
    # the heaviest entry (top) is pinned to the transmuted scale ~173 GeV
    assert np.isclose(s["rescaled_top_GeV"], 173.0, rtol=0.05)


def test_rescaled_matrix_is_no_longer_planck():
    mat, M_Pl, Lam = sg.rescaled_matrix_GeV()
    # every entry drops below 1 TeV -- the scale gap is closed
    for t in mat:
        assert np.all(mat[t] < 1e3)


def test_hierarchy_is_exp_of_order_forty():
    s = sg.scale_summary()
    assert 35 < s["ln_hierarchy"] < 42      # exp(-~40) ~ 10^-17
