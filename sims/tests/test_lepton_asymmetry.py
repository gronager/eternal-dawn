"""Tests for the relic lepton-asymmetry ("unmatched neutrino number") prediction."""

import math

from cartasis_sims import lepton_asymmetry as la


def test_sphaleron_coefficients_are_the_standard_model_values():
    c_B, c_L = la.sphaleron_BL_coefficients()
    assert math.isclose(c_B, 28.0 / 79.0, rel_tol=1e-12)
    assert math.isclose(c_L, -51.0 / 79.0, rel_tol=1e-12)


def test_lepton_asymmetry_is_order_baryon_asymmetry_not_enhanced():
    # sphalerons leave |L/B| = 51/28 ~ 1.82 -- comparable, NOT a large enhancement
    assert math.isclose(la.lepton_to_baryon_ratio(), 51.0 / 28.0, rel_tol=1e-12)
    assert 1.0 < la.lepton_to_baryon_ratio() < 3.0


def test_predicted_degeneracy_is_of_order_eta():
    # xi_nu ~ eta ~ 1e-9, indistinguishable from zero
    xi = la.predicted_neutrino_degeneracy()
    assert 1e-10 < xi < 1e-8


def test_prediction_sits_far_below_the_observational_ceiling():
    assert la.predicted_neutrino_degeneracy() < la.XI_CEILING
    # ~1e7 of unprobed headroom that SCT predicts is empty
    assert la.unmatched_headroom() > 1e5
    assert la.is_discriminating()


def test_degeneracy_inversion_is_linear_in_small_asymmetry():
    # doubling the lepton asymmetry doubles xi in the small-xi regime
    x1 = la.degeneracy_from_lepton_asymmetry(1e-9)
    x2 = la.degeneracy_from_lepton_asymmetry(2e-9)
    assert math.isclose(x2 / x1, 2.0, rel_tol=1e-9)
