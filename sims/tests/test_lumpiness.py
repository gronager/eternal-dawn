"""Tests for the lumpiness ladder across Cartasis membranes."""

import math

import numpy as np

from cartasis_sims import lumpiness as lp


def test_ladder_is_monotonic_increasing():
    L = lp.lumpiness_ladder()
    assert lp.is_monotonic(L)                     # more lumpiness deeper in
    assert L[0] < L[1] < L[2]                      # void-ish OGU < BHU1 < BHU2


def test_biggest_jump_is_ogu_to_bhu1():
    # the parent's lumpiness first appears at the OGU -> BHU1 crossing
    L = lp.lumpiness_ladder()
    assert lp.biggest_jump_index(L) == 0


def test_ladder_converges_to_fixed_point():
    a, b = 1.0, 0.7
    L = lp.lumpiness_ladder(a=a, b=b, n_max=30)
    Lstar = lp.fixed_point(a, b)
    assert math.isclose(L[-1], Lstar, rel_tol=1e-3)
    assert math.isclose(Lstar, a / math.sqrt(1 - b**2), rel_tol=1e-9)


def test_ogu_has_minimal_lumpiness():
    # the OGU (LambdaCDM smooth-start analogue) is the low rung, well below the BHUs
    L = lp.lumpiness_ladder()
    assert L[0] < 0.5 * L[1]


def test_ultra_large_structures_exceed_homogeneity():
    # the named ultra-large structures are all above the homogeneity scale; BAO is not
    for name, mpc in lp.ULTRA_LARGE.items():
        if name == "BAO scale (ref)":
            assert not lp.exceeds_homogeneity(mpc)
        else:
            assert lp.exceeds_homogeneity(mpc)


def test_gly_to_k_conversion():
    # 1 Gly ~ 307 Mpc -> k ~ 0.02 /Mpc
    assert math.isclose(lp.gly_to_k(1.0), 2 * math.pi / 306.6, rel_tol=1e-6)


def test_lcdm_floor_is_tiny_on_ultra_large_scales():
    # CAMB-computed: sigma8 ~ 0.81 but the floor collapses to ~0.01 at the
    # homogeneity scale -- so observed ultra-large structure is genuine excess.
    import pytest
    if not lp.HAVE_CAMB:
        pytest.skip("CAMB not installed")
    assert abs(lp.lcdm_sigma_R(8.0) - 0.81) < 0.05      # sigma8 anchor
    assert lp.lcdm_sigma_R(260.0) < 0.05                 # floor tiny at homogeneity
    # monotonically falling with scale
    assert lp.lcdm_sigma_R(400.0) < lp.lcdm_sigma_R(100.0)
