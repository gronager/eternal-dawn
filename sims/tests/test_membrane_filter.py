"""Tests for the membrane-filter dark-to-baryon ratio f ~ 1/6.

Includes a SymPy check of the pass-fraction inversions -- the algebra the whole
'f from an O(1) torsion barrier' argument rests on.
"""

import math

import pytest

from cartasis_sims import membrane_filter as mf


def test_observed_f_is_about_one_sixth():
    f = mf.observed_f(0.02237, 0.1200)
    assert abs(f - 1.0 / 6.0) < 0.02
    assert abs(f - 0.157) < 0.005


def test_pass_fraction_inversions_round_trip():
    for kind in ("boltzmann", "fermi"):
        for x in (0.5, 1.0, 1.79, 2.5):
            f = (mf.pass_fraction_boltzmann(x) if kind == "boltzmann"
                 else mf.pass_fraction_fermi(x))
            assert math.isclose(mf.barrier_ratio_for_f(f, kind), x, rel_tol=1e-9)


def test_one_sixth_maps_to_ln6_and_ln5():
    # f = 1/6 <=> x = ln6 (Boltzmann), x = ln5 (Fermi)
    assert math.isclose(mf.barrier_ratio_for_f(1 / 6, "boltzmann"), math.log(6), rel_tol=1e-9)
    assert math.isclose(mf.barrier_ratio_for_f(1 / 6, "fermi"), math.log(5), rel_tol=1e-9)


def test_observed_f_needs_order_unity_barrier():
    # the measured f sits at x ~ 1.7-1.9 -- an O(1) barrier, not 1e-9 or 1e9
    f = mf.observed_f()
    xB = mf.barrier_ratio_for_f(f, "boltzmann")
    xF = mf.barrier_ratio_for_f(f, "fermi")
    assert 1.5 < xB < 2.0
    assert 1.4 < xF < 1.9


def test_one_sixth_is_natural_for_an_order_unity_barrier():
    # an O(1) torsion barrier (x in [1, 2.5]) produces f in the ~0.08-0.37 decade,
    # bracketing the observed 1/6: it is natural, not fine-tuned
    for kind in ("boltzmann", "fermi"):
        f_min, f_max = mf.natural_band(kind)
        assert f_min < mf.observed_f() < f_max
        assert mf.is_f_natural(kind=kind)


def test_sympy_confirms_the_inversions():
    sp = pytest.importorskip("sympy")
    f, x = sp.symbols("f x", positive=True)
    xB = sp.solve(sp.exp(-x) - f, x)[0]
    xF = sp.solve(1 / (sp.exp(x) + 1) - f, x)[0]
    # symbolic forms
    assert sp.simplify(xB - (-sp.log(f))) == 0
    assert sp.simplify(xF - sp.log((1 - f) / f)) == 0
    # f = 1/6 special values, symbolically
    assert sp.simplify(xB.subs(f, sp.Rational(1, 6)) - sp.log(6)) == 0
    assert sp.simplify(xF.subs(f, sp.Rational(1, 6)) - sp.log(5)) == 0
