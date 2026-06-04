"""Tests for the Weltformel physical-units spectrum (all 15 masses, no fitting)."""

import numpy as np

from cartasis_sims import weltformel as wf


def test_structure_is_exact():
    mat = wf.fermion_matrix()
    # gen III heaviest in every tower; quarks above leptons at every generation
    for t in wf.TOWERS:
        assert mat[t][2] >= mat[t][1] >= mat[t][0]
    for i in range(3):
        assert mat["up-quark"][i] > mat["charged-lepton"][i]
        assert mat["charged-lepton"][i] >= mat["neutrino"][i]


def test_neutrino_leading_dirac_vanishes():
    # the right-handed neutrino is a total singlet, so the LEADING Dirac mass vanishes;
    # the physical mass is the (tiny, nonzero) sub-eV seesaw floor, with the scale M owed
    mat = wf.fermion_matrix()
    assert np.allclose(mat["neutrino"], 0.0)        # leading order
    # and the leading row is far below the charged leptons (seesaw-suppressed, not comparable)
    assert mat["neutrino"][2] < 1e-6 * mat["charged-lepton"][2]


def test_single_global_scale_pins_the_top():
    # the one anchor is the heaviest fermion = the condensate scale
    mat = wf.fermion_matrix(top_GeV=173.0)
    assert np.isclose(mat["up-quark"][2], 173.0, rtol=1e-6)


def test_charged_fermions_within_two_orders():
    # no per-tower fitting, so the worst case is ~1-2 orders -- but every charged fermion is
    # placed within that, far better than 'unknown'
    assert wf.summary()["worst_fermion_factor"] < 100


def test_bosons_within_a_few_percent():
    assert wf.summary()["worst_boson_factor"] < 1.05
