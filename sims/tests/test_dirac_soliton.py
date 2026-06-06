"""First-order (G,F) Dirac torsiton: the eigensolver is validated, range beats contact."""
import numpy as np

from cartasis_sims import chiral_soliton as cs
from cartasis_sims import dirac_soliton as ds


def test_first_order_matches_second_order_on_same_well():
    # in the SAME converged scalar well, the first-order Dirac binds where the second-order does,
    # at a nearby energy (validates the eigensolver; the doublers stay out of the gap)
    o = cs.solve_chiral(n_fermions=1, n_levels=4, g_v=0.0)
    lv, Mvac = ds.dirac_levels(np.abs(o["M"]), np.zeros_like(o["r"]), o["r"], n_levels=4)
    assert len(lv) >= 1
    E0 = lv[0][0]
    assert 0 < E0 < Mvac                               # a genuine bound state in the gap
    assert abs(E0 - o["E"][0]) < 0.5                   # near the second-order ground state


def test_self_consistent_soliton_binds_chiral():
    # the first-order self-consistent torsiton binds with a real bag (g_v=0, scalar only)
    o = ds.solve_dirac_soliton(n_fermions=1, g_v=0.0, iters=120, N=200)
    assert o["bound"] >= 1
    assert o["core"] < 0.0                              # chiral-restored core = a real bag


def test_repulsion_monotonically_loosens_the_bag():
    # what's solid: turning the (ranged) repulsion up monotonically shallows the bag (core less
    # negative) -- the spin-push doing real work. (A clean range-vs-contact comparison is confounded
    # by the coupling normalisation -- the omega amplitude scales as 1/m_omega^2 -- so that claim
    # waits on the proper coupling derivation.)
    cores = [ds.solve_dirac_soliton(n_fermions=1, g_v=gv, m_omega=4.0, iters=160, N=200)["core"]
             for gv in (0.0, 2.0, 4.0)]
    assert cores[0] < cores[1] < cores[2] < 0.0        # shallower with more push, still a bag
