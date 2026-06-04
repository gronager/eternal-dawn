"""Tests for the 1-D flux-tube (Nielsen-Olesen vortex) string tension."""

import math

import numpy as np

from cartasis_sims import flux_tube as ft


def test_vortex_converges_with_right_boundary_conditions():
    out = ft.solve_vortex(beta=2.0)
    assert out["success"]
    # condensate: 0 on the axis (normal core), -> 1 outside (superconductor)
    assert abs(out["f"][0]) < 0.05
    assert abs(out["f"][-1] - 1.0) < 0.02
    # gauge profile: 0 on axis -> 1 outside (flux quantised)
    assert abs(out["a"][0]) < 0.05
    assert abs(out["a"][-1] - 1.0) < 0.02


def test_bps_tension_is_2pi():
    # topological saturation at the BPS point (beta=2, m_H=m_W): sigma = 2 pi v^2
    assert math.isclose(ft.bps_tension(), 2 * math.pi, rel_tol=2e-3)


def test_flux_is_localized_in_a_tube():
    out = ft.solve_vortex(beta=2.0)
    B = np.abs(out["B"])
    # the colour-electric field is concentrated near the axis and falls off
    assert B[5] > B[len(B) // 2] > B[-1]


def test_tension_rises_with_beta():
    # type-I (beta<2) below the BPS value, type-II (beta>2) above
    assert ft.string_tension(1.0) < ft.bps_tension() < ft.string_tension(3.0)


def test_confining_potential_is_linear():
    L = np.array([1.0, 2.0, 4.0])
    V = ft.confining_potential(L, sigma=2 * math.pi)
    assert np.allclose(V / L, 2 * math.pi)        # V proportional to L (confinement)
