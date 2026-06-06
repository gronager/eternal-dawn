"""The Fierz vector/axial repulsion in the chiral soliton -- it matters, and its naive coupling
unbinds the soliton (a normalisation crux, documented honestly)."""
import numpy as np

from cartasis_sims import chiral_soliton as cs


def _bound(o):
    return int(np.sum(o["E"] < o["constituent_mass"]))


def test_default_is_scalar_only_and_binds():
    # g_v=0 default preserves the prior scalar-only soliton: it binds
    o = cs.solve_chiral(n_fermions=3, n_levels=6, iters=600, mix=0.15, g_v=0.0)
    assert o["converged"] and _bound(o) >= 1


def test_repulsion_reduces_binding_monotonically():
    # turning the vector repulsion up (below the unbinding threshold) monotonically loosens the bag
    masses = [cs.solve_chiral(n_fermions=3, n_levels=6, iters=600, mix=0.15, g_v=gv)["mass"]
              for gv in (0.0, 1.0, 2.0, 3.0)]
    assert masses[0] > masses[1] > masses[2] > masses[3]   # the repulsion does real work


def test_fierz_naive_coupling_unbinds_the_soliton():
    # THE FINDING: g_v = g*sqrt(2) (the Fierz-naive value) overwhelms the attraction -> no bound
    # state. Survival needs g_v/g <~ 1; the naive mapping is a normalisation artefact, the true
    # coupling is an open derivation. This test pins the behaviour so it isn't silently 'fixed'.
    g = 4.0
    o_survive = cs.solve_chiral(n_fermions=3, n_levels=6, iters=600, mix=0.15, g_v=0.75 * g)
    o_unbind = cs.solve_chiral(n_fermions=3, n_levels=6, iters=600, mix=0.15, g_v=cs.fierz_gv(g))
    assert _bound(o_survive) >= 1                          # below threshold: soliton survives
    assert _bound(o_unbind) == 0                           # Fierz-naive: unbound
    assert np.isclose(cs.fierz_gv(g), g * np.sqrt(2.0))
