"""Tests for the Phase-0 black-hole identities."""

import math

from cartasis_sims import blackhole as bh
from cartasis_sims import constants as k


def test_schwarzschild_equals_hubble_when_flat():
    # R_s / R_H = Omega; for a flat universe it is exactly 1.
    assert math.isclose(bh.rs_over_rh(Omega=1.0), 1.0, rel_tol=1e-12)
    assert math.isclose(bh.rs_over_rh(Omega=0.7), 0.7, rel_tol=1e-12)


def test_max_spin_is_order_two():
    # a*_max = 2 / Omega.
    assert math.isclose(bh.max_spin_parameter(Omega=1.0), 2.0, rel_tol=1e-12)


def test_hubble_mass_order_1e53():
    M = bh.mass_within_hubble(Omega=1.0)
    assert 5e52 < M < 2e53


def test_hawking_is_thirty_orders_below_cmb():
    M = bh.mass_within_hubble(Omega=1.0)
    T_H = bh.hawking_temperature(M)
    assert 1e-31 < T_H < 1e-29
    assert k.T_cmb / T_H > 1e29


def test_hawking_equals_half_gibbons_hawking():
    # Identity that holds when R_s = R_H.
    M = bh.mass_within_hubble(Omega=1.0)
    assert math.isclose(bh.hawking_temperature(M),
                        bh.gibbons_hawking_temperature() / 2.0, rel_tol=1e-9)


def test_filter_fraction_one_sixth():
    f = bh.filter_fraction()
    assert math.isclose(f, k.Omega_b_hsq / (k.Omega_b_hsq + k.Omega_c_hsq))
    assert 0.14 < f < 0.18
