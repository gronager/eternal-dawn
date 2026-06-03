"""Tests for the colour-channel forces and the residual overlap force."""

import math

import numpy as np

from cartasis_sims import color_force as cf


def test_fundamental_casimir():
    # Gell-Mann validation: <T.T> = C_F = 4/3 for the fundamental
    assert math.isclose(cf.casimir_fundamental(), 4.0 / 3.0, rel_tol=1e-9)


def test_qcd_color_factors():
    assert math.isclose(cf.color_factor("qqbar_singlet"), -4.0 / 3.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qqbar_octet"), 1.0 / 6.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qq_antitriplet"), -2.0 / 3.0, abs_tol=1e-9)
    assert math.isclose(cf.color_factor("qq_sextet"), 1.0 / 3.0, abs_tol=1e-9)


def test_only_antisymmetric_channels_bind():
    # singlet and antitriplet attract; octet and sextet repel
    assert cf.binds("qqbar_singlet")
    assert cf.binds("qq_antitriplet")
    assert not cf.binds("qqbar_octet")
    assert not cf.binds("qq_sextet")


def test_singlet_is_most_attractive():
    # the meson (q-qbar singlet) is more attractive than the diquark (antitriplet)
    assert cf.color_factor("qqbar_singlet") < cf.color_factor("qq_antitriplet") < 0


def test_residual_force_is_attractive_and_short_range():
    r = np.array([0.5, 1.0, 2.0, 4.0])
    V = cf.residual_yukawa(r)
    assert np.all(V < 0)                       # attractive
    assert abs(V[-1]) < abs(V[0])              # falls off (short range)
