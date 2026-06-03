"""Tests for the tensor-to-scalar ratio / inflation-vs-bounce discriminator."""

import math

from cartasis_sims import tensor_ratio as tr


def test_starobinsky_lands_at_planck_ns_low_r():
    n_s, r = tr.starobinsky(57.0)
    assert abs(n_s - 0.965) < 0.003          # plateau n_s ~ 0.965
    assert 0.002 < r < 0.006                  # r ~ 0.0037, within LiteBIRD reach


def test_chaotic_phi2_is_excluded():
    n_s, r = tr.chaotic_phi2(57.0)
    assert r > tr.BK18_UPPER                  # r ~ 0.14 > 0.036, ruled out


def test_inflation_consistency_relation():
    # n_T = -r/8 for slow-roll
    assert math.isclose(tr.inflation_consistency_nT(0.04), -0.005)


def test_bounce_equal_tilt_relation():
    # the matter bounce predicts n_T = n_s - 1 (tilts coincide)
    assert math.isclose(tr.bounce_tensor_tilt(0.9649), -0.0351, abs_tol=1e-4)


def test_bounce_tilt_from_w_matches_equal_tilt():
    # the w that gives n_s=0.965 also gives n_T = n_s-1 via the MS evolution
    w = tr.w_for_observed_tilt(0.9649)
    nT = tr.bounce_tensor_tilt_from_w(w)
    assert math.isclose(nT, tr.bounce_tensor_tilt(0.9649), abs_tol=2e-3)


def test_discriminator_gap_is_large():
    # at the observed n_s the two hypotheses' tensor tilts differ by ~70x
    d = tr.discriminator_gap()
    assert d["nT_bounce"] < -0.03            # bounce: clearly red
    assert abs(d["nT_inflation"]) < 0.001    # plateau inflation: ~flat
    assert abs(d["gap"]) > 0.03              # a wide, in-principle-measurable gap
