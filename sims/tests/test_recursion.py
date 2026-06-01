"""Tests for the OGU-as-finite-recursive-timeline picture."""

import math

from cartasis_sims import recursion as rc


def test_equilibrium_mass_is_a_hubble_mass():
    # M_eq where T_BH = T_dS is ~ a Hubble mass (~ our universe), dark-energy set
    M = rc.equilibrium_mass()
    assert 1e52 < M < 1e53
    # at M_eq the Hawking and de Sitter temperatures match
    assert math.isclose(rc.hawking_temperature(M), rc.de_sitter_temperature(),
                        rel_tol=1e-9)


def test_growth_caps_at_nariai():
    M_N = rc.nariai_mass()
    assert 1e52 < M_N < 1e53                   # ~ a Hubble mass (geometric cap)
    assert not rc.saturates(0.5 * M_N)         # still growing
    assert rc.saturates(1.5 * M_N)             # past the geometric cap


def test_membrane_lives_astronomically_long_but_finite():
    t = rc.membrane_lifetime()
    assert t > 1e140                            # ~1e142 s -- not eternal, but vast
    assert math.isfinite(t)


def test_many_recursion_cycles_per_lifetime():
    # the interior runs ~1e34 full recursive cycles before the membrane evaporates
    n = rc.recursion_cycles_per_lifetime()
    assert 1e32 < n < 1e36


def test_never_touches_a_neighbour():
    # dilute regime: neighbours are e^{I/4} horizons away, growth caps at ~1 horizon
    assert not rc.touches_neighbour_before_evaporating(1e80)
    assert rc.touches_neighbour_before_evaporating(0.5)   # only if essentially adjacent


def test_hawking_temperature_falls_with_mass():
    assert rc.hawking_temperature(1e52) > rc.hawking_temperature(1e53)
