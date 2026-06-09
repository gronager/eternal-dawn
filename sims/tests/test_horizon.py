"""The Horizon spectrum: one condensate, all twelve fermions, the hierarchy factorising as c_T x rung."""
import numpy as np
from cartasis_sims import horizon as hz


def test_relative_masses_factorise():
    rel = hz.relative_masses()
    assert set(rel) == set(hz.TOWERS)
    # charged leptons relative to the electron: e=1, mu~207, tau~3477
    m = rel["charged-lepton"]["m"]
    assert np.isclose(m[0], 1.0) and abs(m[1] - 206.8) < 1 and abs(m[2] - 3477) < 5
    # the framework's claim: m = c_T * rung exactly (factorisation is consistent)
    for t in hz.TOWERS:
        assert np.allclose(rel[t]["c_T"] * rel[t]["rungs"], rel[t]["m"], rtol=1e-9)
    # colour/charge ordering: up-type heaviest tower, neutrino lightest
    cT = {t: rel[t]["c_T"] for t in hz.TOWERS}
    assert cT["up-quark"] > cT["down-quark"] > cT["charged-lepton"] > cT["neutrino"]


def test_emergence_switches_on_at_one_transition():
    temps, v, masses, rel = hz.emergence(n_T=200, Tc=0.5, width=0.05)
    assert v[0] < 0.02 and v[-1] > 0.98                  # condensate off (hot) -> on (cold)
    # every tower is ~massless hot and at its relative mass cold
    for t in ("up-quark", "down-quark", "charged-lepton"):
        assert masses[t][0].max() < 0.05 * rel[t]["m"].max()    # hot: suppressed
        assert np.allclose(masses[t][-1], rel[t]["m"], rtol=1e-6)  # cold: full hierarchy
