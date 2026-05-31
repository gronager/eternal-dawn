"""Tests for the axis-alignment machinery (Corner B)."""

import numpy as np

from cartasis_sims import axes as ax


def test_lb_to_vec_pole_and_origin():
    assert np.allclose(ax.lb_to_vec(0, 90), [0, 0, 1], atol=1e-12)
    assert np.allclose(ax.lb_to_vec(0, 0), [1, 0, 0], atol=1e-12)
    assert np.allclose(ax.lb_to_vec(90, 0), [0, 1, 0], atol=1e-12)


def test_acute_angle_is_sign_independent():
    v = ax.lb_to_vec(123, 45)
    assert ax.acute_angle_deg(v, v) < 1e-3             # arccos near 1: sqrt-precision
    assert ax.acute_angle_deg(v, -v) < 1e-3            # antipode = same axis
    assert abs(ax.acute_angle_deg(ax.lb_to_vec(0, 0),
                                  ax.lb_to_vec(90, 0)) - 90.0) < 1e-9


def test_axis_pvalue_endpoints():
    assert ax.axis_pvalue(0.0) == 0.0
    assert abs(ax.axis_pvalue(90.0) - 1.0) < 1e-12
    assert abs(ax.axis_pvalue(60.0) - 0.5) < 1e-12     # 1 - cos60 = 0.5


def test_concentration_aligned_vs_isotropic():
    aligned = np.array([ax.lb_to_vec(10, 80), ax.lb_to_vec(190, 82),
                        ax.lb_to_vec(30, 78)])
    tau1_a, _ = ax.concentration(aligned)
    # three mutually orthogonal axes are maximally isotropic -> tau1 = 1/3
    iso = np.eye(3)
    tau1_i, _ = ax.concentration(iso)
    assert tau1_a > 0.9
    assert abs(tau1_i - 1.0 / 3.0) < 1e-9


def test_monte_carlo_calibration():
    # Tightly clustered axes must be improbable under isotropy; an orthogonal
    # (maximally spread) triad must not be.
    clustered = np.array([ax.lb_to_vec(10, 80), ax.lb_to_vec(350, 84),
                          ax.lb_to_vec(40, 86)])
    p_clustered = ax.monte_carlo_pvalue(clustered, n_trials=20_000, seed=3)["p_value"]
    p_ortho = ax.monte_carlo_pvalue(np.eye(3), n_trials=20_000, seed=3)["p_value"]
    assert p_clustered < 0.05
    assert p_ortho > 0.5
