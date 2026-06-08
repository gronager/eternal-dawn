"""The Woods-Saxon bag model of the torsiton tower: exactly 3 rungs, 0/1/2 nodes, the cap."""
import numpy as np
from cartasis_sims import woods_saxon as ws


def test_three_rungs_with_hydrogen_node_structure():
    r, states = ws.ws_spectrum(V0=5.0, R0=3.0, a=0.5)
    assert len(states) == 3                                   # exactly three generations bind
    E = [e for e, _ in states]
    assert E[0] < E[1] < E[2] < 0                             # bound, ordered
    assert [ws.interior_nodes(u) for _, u in states] == [0, 1, 2]   # 0/1/2 nodes


def test_a_fourth_rung_appears_only_when_deeper():
    # the cap is real: deepening the well past threshold binds a fourth (so V0=5 is "three, not four")
    assert len(ws.ws_spectrum(V0=5.0, R0=3.0, a=0.5)[1]) == 3
    assert len(ws.ws_spectrum(V0=7.0, R0=3.0, a=0.5)[1]) == 4


def test_overlap_mass_rises_with_generation():
    # Eq. configmass weights the density by the LOCAL mass M(r), large in the VACUUM (M = 1 - sigma).
    # The ground rung sits in the core where M is small (lightest); higher rungs spread into the
    # vacuum where M is large -> heavier. This is the physical ordering (electron < mu < tau).
    r, states = ws.ws_spectrum(V0=5.0, R0=3.0, a=0.5)
    cond = 1.0 / (1.0 + np.exp((r - 3.0) / 0.5))              # condensate sigma(r), large inside the bag
    mass_r = 1.0 - cond                                       # local mass M(r), large in the vacuum
    m = ws.overlap_masses(r, states, mass_r)
    assert 0 < m[0] < m[1] < m[2]
    # and weighting by sigma instead would give the WRONG (falling) order -- the catch that fixed this
    msig = ws.overlap_masses(r, states, cond)
    assert msig[0] > msig[1] > msig[2] > 0
