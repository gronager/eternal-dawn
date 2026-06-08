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


def test_overlap_mass_falls_with_generation():
    # mass = overlap with the condensate; higher rungs spread/nodal -> less overlap -> lighter
    r, states = ws.ws_spectrum(V0=5.0, R0=3.0, a=0.5)
    cond = 1.0 / (1.0 + np.exp((r - 3.0) / 0.5))
    m = ws.overlap_masses(r, states, cond)
    assert m[0] > m[1] > m[2] > 0
