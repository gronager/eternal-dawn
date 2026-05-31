"""Tests for the recursive cycle (Weyl entropy reset, de Sitter void)."""

import math

from cartasis_sims import cycle as cy


def test_weyl_entropy_cycles_low_high_low():
    # gravitational (clumping) entropy: low at bounce, high at black-hole peak,
    # low again at de Sitter heat death -- the reset that allows re-nucleation
    s_bounce = cy.weyl_entropy(cy.T_BOUNCE)
    s_peak = cy.weyl_entropy(3e40)            # mid-life, holes present
    s_dead = cy.weyl_entropy(cy.T_DESITTER)
    assert s_peak > s_bounce + 5.0
    assert s_peak > s_dead + 5.0
    assert abs(s_dead - s_bounce) < 1.0       # resets to ~ its bounce value


def test_total_entropy_is_monotonic():
    import numpy as np
    t = np.logspace(-36, 112, 400)
    S = cy.total_entropy(t)
    assert np.all(np.diff(S) >= -1e-9)        # Second Law: never decreases
    assert S[-1] > S[0] + 25.0                # rises to the de Sitter maximum


def test_de_sitter_temperature_is_about_1e_minus_30_K():
    T = cy.de_sitter_temperature(1.8e-18)     # our Lambda
    assert 1e-31 < T < 1e-29


def test_cycle_always_closes_given_infinite_time():
    assert cy.cycle_closes(0.0)               # even a zero-ish rate fires eventually
    assert cy.cycle_closes(1e-300)
