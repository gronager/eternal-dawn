"""Tests for the Tier-1 Einstein-Cartan bounce."""

import numpy as np

from cartasis_sims import bounce as bnc


def test_bounce_reaches_rho_C_and_no_singularity():
    sol = bnc.simulate_bounce(w=1.0 / 3.0, a_init=3.0)
    # bounce at a_min = 1 where rho = rho_C
    assert abs(sol.a_min - 1.0) < 1e-3
    assert abs(sol.rho_max - 1.0) < 1e-3
    # never reaches a singularity
    assert sol.a.min() > 0.5


def test_hubble_changes_sign_at_bounce():
    sol = bnc.simulate_bounce(w=1.0 / 3.0, a_init=3.0)
    i = int(np.argmin(sol.a))
    assert sol.H[0] < 0          # starts collapsing
    assert sol.H[-1] > 0         # ends expanding
    assert abs(sol.H[i]) < 1e-3  # H ~ 0 at the bounce


def test_bounce_is_time_symmetric():
    sol = bnc.simulate_bounce(w=1.0 / 3.0, a_init=3.0)
    # density profile should be symmetric about the bounce (t -> -t)
    i = int(np.argmin(sol.a))
    span = min(i, len(sol.a) - 1 - i)
    left = sol.rho[i - span:i]
    right = sol.rho[i + 1:i + 1 + span][::-1]
    assert np.allclose(left, right, rtol=1e-2, atol=1e-2)


def test_gr_collapse_hits_zero():
    gr = bnc.simulate_gr_collapse(w=1.0 / 3.0, a_init=3.0)
    assert gr["t_singularity"] > 0
    assert gr["a"][-1] < 0.1 * 3.0   # dives toward a = 0 (well below a_init)


def test_physical_timescale_is_femtoscale():
    tau = bnc.physical_timescale(1.0e50)
    assert 1e-22 < tau < 1e-19   # ~10^-21 s for rho_C ~ 1e50 kg/m^3
