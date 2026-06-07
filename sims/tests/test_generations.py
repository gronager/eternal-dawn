"""Tests for the generations overlap-hierarchy mechanism."""

import numpy as np

from cartasis_sims import generations as gen
from cartasis_sims import soliton as so


def test_sm_fermion_masses_are_approximately_geometric():
    # the observed fermion masses are evenly spaced in log to within a factor ~2
    for name, m in gen.SM_FERMIONS.items():
        assert gen.is_approximately_geometric(m), name
        gaps = gen.log_spacings(m)
        assert np.all(gaps > 0)                       # monotonic hierarchy


def test_overlap_falls_with_level():
    lv = so.energy_levels(n_levels=4, kind="linear", depth=6.0, E_max=18.0, n_scan=100)
    O = np.sqrt(gen.overlap_masses(lv, source_size=0.25))
    assert np.all(np.diff(O) < 0)                     # higher levels overlap the core less


def test_overlap_hierarchy_is_geometric():
    # ln(mass) should be ~linear in level index (a geometric hierarchy)
    lv = so.energy_levels(n_levels=4, kind="linear", depth=6.0, E_max=18.0, n_scan=100)
    m = gen.overlap_masses(lv, source_size=0.3)
    n = np.arange(len(m))
    resid = np.log(m) - np.polyval(np.polyfit(n, np.log(m), 1), n)
    rms = np.sqrt(np.mean(resid**2)) / abs(np.mean(np.log(m)))
    assert rms < 0.15                                 # log-linear (geometric) to ~15%


def test_self_consistent_hierarchy_broad_vs_sharp():
    # the real (self-consistent-well) calculation: arithmetic levels, broad overlap = NO hierarchy,
    # sharp overlap = a geometric hierarchy, with the light generation anomalously light
    h = gen.self_consistent_hierarchy(g=14.0, m_sigma=0.5)
    assert h["converged"] and h["n_bound"] >= 3
    # level energies roughly arithmetic (spacings within a factor ~2)
    s = np.abs(h["level_spacings"])
    assert s.max() / s.min() < 2.5
    # broad substrate gives essentially no hierarchy; the sharp core gives a real one
    assert h["broad_spread"] < 2.0
    assert h["sharp_spread"] > 10.0
    m = h["sharp_masses"]
    assert m[0] < m[1] < m[2]                              # monotonic, geometric tower
    # the big jump is at the light end (the "electron is anomalously light" shape), as observed
    assert np.log(m[1] / m[0]) > np.log(m[2] / m[1])


def test_steeper_source_gives_bigger_factor():
    # tightening the core source localises the mass-source -> steeper hierarchy
    lv = so.energy_levels(n_levels=4, kind="linear", depth=6.0, E_max=18.0, n_scan=100)
    f_small = gen.geometric_factor(gen.overlap_masses(lv, source_size=0.15))
    f_big = gen.geometric_factor(gen.overlap_masses(lv, source_size=0.4))
    assert f_big > f_small > 1.0
