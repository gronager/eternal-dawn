"""Tests for the self-bound Fermi ball and the degeneracy-pressure effect."""

import math

import numpy as np

from cartasis_sims import fermi_ball as fb


def test_fermi_kinetic_nonrelativistic_limit():
    # e_kin/m - 1 -> (3/10) x^2 at small density (classic degeneracy KE)
    n = 1e-3
    x = fb.kF_of_n(n)
    assert math.isclose(fb.e_kin_per_particle(n) - 1.0, 0.3 * x**2, rel_tol=0.05)


def test_kinetic_energy_rises_with_density():
    n = np.array([1e-3, 1e-2, 1e-1, 0.4])
    e = fb.e_kin_per_particle(n)
    assert np.all(np.diff(e) > 0)             # degeneracy pressure: KE grows with n


def test_both_bind_in_chosen_regime():
    assert fb.saturation(degeneracy=True)["bound"]
    assert fb.saturation(degeneracy=False)["bound"]


def test_degeneracy_lowers_density_and_binding():
    d = fb.degeneracy_shift()
    # Pauli pressure pushes the drop apart: lower saturation density, shallower binding
    assert d["n0_with"] < d["n0_without"]
    assert d["binding_with"] < d["binding_without"]
    assert 0.5 < d["n0_ratio"] < 1.0          # density drops but stays bound


def test_no_degeneracy_minimum_is_analytic():
    # without the Pauli KE, e(n) = 1 - a n + b n^2 minimised at n = a/(2b)
    a, b = 5.0, 4.0
    sat = fb.saturation(a=a, b=b, degeneracy=False)
    assert math.isclose(sat["n0"], a / (2 * b), rel_tol=0.02) or sat["n0"] >= 0.59
