"""Tests for the ab-initio spectrum from G, hbar, c alone (the Weltformel computation)."""

import numpy as np

from cartasis_sims import ab_initio_spectrum as sp


def test_planck_mass_from_constants():
    # M_Pl c^2 ~ 1.22e19 GeV from G, hbar, c -- the only scale
    assert np.isclose(sp.planck_mass_GeV(), 1.22e19, rtol=2e-2)


def test_matrix_is_3x4_and_planck_scaled():
    fm = sp.fermion_matrix_GeV()
    assert set(fm["matrix_GeV"]) == set(sp.TOWERS)
    for t in sp.TOWERS:
        row = fm["matrix_GeV"][t]
        assert len(row) == 3
        # every entry is near the Planck scale (10^17-10^19 GeV) -- off by the hierarchy
        assert np.all(row > 1e16) and np.all(row < 1e20)


def test_structure_ordering_quarks_above_leptons_above_neutrino():
    fm = sp.fermion_matrix_GeV()
    # at fixed generation, quarks > charged lepton > neutrino (charge/colour handles, no fit)
    for i in range(3):
        assert fm["matrix_GeV"]["up-quark"][i] > fm["matrix_GeV"]["charged-lepton"][i]
        assert fm["matrix_GeV"]["charged-lepton"][i] > fm["matrix_GeV"]["neutrino"][i]


def test_generation_ordering_heaviest_is_gen_three():
    fm = sp.fermion_matrix_GeV()
    for t in sp.TOWERS:
        row = fm["matrix_GeV"][t]
        assert row[2] > row[1] and row[1] > row[0]   # gen III heaviest (ground state)


def test_hierarchy_gap_is_enormous():
    # the one number not produced: the Planck/electroweak ratio ~10^16-17
    assert sp.hierarchy_gap() > 1e15


def test_bosons_are_also_planck_scaled():
    bm = sp.boson_masses_GeV()
    for k in ("W", "Z", "H"):
        assert bm[k] > 1e15            # Planck-scaled like the fermions
    assert bm["Z"] > bm["W"]           # custodial ordering preserved
