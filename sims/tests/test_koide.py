"""Koide: charged leptons hit 2/3 to 1e-4; the Z3 fit puts the electron at the overlap near-node."""
import numpy as np

from cartasis_sims import koide as kd


def test_charged_leptons_two_thirds():
    Q = kd.koide_ratio(kd.LEPTONS)
    assert abs(Q - 2 / 3) < 1e-4               # the famous 2/3 to 1e-5


def test_z3_fit_reproduces_lepton_masses():
    m0, d, perm, pred = kd.z3_fit(kd.LEPTONS)
    # two parameters (m0, delta), sqrt2 amplitude fixed, reproduce all three to a few percent
    assert np.allclose(np.sort(pred), np.sort(kd.LEPTONS), rtol=0.05)


def test_electron_is_the_overlap_near_node():
    # the hierarchy lives in the phase: the lightest rung sits where (1+sqrt2 cos) ~ 0 (overlap cancels)
    m0, d, perm, pred = kd.z3_fit(kd.LEPTONS)
    fac = kd.node_proximity(m0, d, perm)
    assert np.min(np.abs(fac)) < 0.1           # one sector is near the node = anomalously light
    # and that near-node sector is the lightest predicted mass
    assert np.argmin(np.abs(fac)) == np.argmin(pred)


def test_quarks_only_approximate():
    # leptons clean (exact Z3), quarks have broken Z3 (CKM mixing) -> Koide only approximate
    assert abs(kd.koide_ratio(kd.UP_QUARKS) - 2 / 3) > 0.1
