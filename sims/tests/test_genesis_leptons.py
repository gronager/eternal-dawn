"""Horizon 3c: the lepton sector -- charged-lepton generations (radial rungs) + suppressed neutrinos."""
import numpy as np
from cartasis_sims import genesis_leptons as lep


def test_charged_lepton_generations():
    c = lep.charged_lepton_tower()
    assert c["names"] == ["e", "mu", "tau"]
    assert abs(c["mass"][2] - c["obs"][2]) < 1.0           # tau anchors the tower
    assert abs(c["mass"][0] - c["obs"][0]) / c["obs"][0] < 0.1   # electron predicted within 10%
    assert c["mass"][0] < c["mass"][1] < c["mass"][2]      # the generation ladder, rising
    assert 0.3 < c["s_T"] < 0.7                            # the bag-sharpness rung parameter


def test_neutrinos_are_the_weakest_grip():
    n = lep.neutrino_tower()
    c = lep.charged_lepton_tower()
    assert n["c_ratio"] < 1e-5                             # far weaker condensate coupling
    assert np.all(n["mass"] < 1e-3 * c["mass"])           # the whole tower sits far below the charged
    s = lep.lepton_spectrum()
    assert "charged" in s and "neutrino" in s
