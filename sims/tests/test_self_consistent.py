"""Tests for the self-consistent (Hartree) gravity-torsion soliton."""

import numpy as np

from cartasis_sims import self_consistent as sc


def test_loop_converges():
    out = sc.solve_soliton(g=4.0, m_sigma=1.0, n_fermions=1, N=400)
    assert out["converged"]
    assert out["residual"] < 1e-6
    # the residual decreases overall (it is a contraction to a fixed point)
    h = out["history"]
    assert h[-1] < h[0]


def test_soliton_is_self_bound():
    out = sc.solve_soliton(g=4.0, m_sigma=1.0, n_fermions=1, N=400)
    assert len(out["bound"]) >= 1            # at least one bound level (E < m0)
    assert np.all(out["bound"] < 1.0)


def test_core_is_chiral_restored():
    # the self-dug well dips the effective mass well below m0 (here through zero)
    out = sc.solve_soliton(g=4.0, m_sigma=1.0, n_fermions=1, N=400)
    assert out["M"][0] < 0.5                  # deep core
    assert out["M"][-1] > 0.9                  # -> m0 outside (localized)


def test_stronger_coupling_binds_deeper():
    weak = sc.solve_soliton(g=3.0, m_sigma=1.0, n_fermions=1, N=400)
    strong = sc.solve_soliton(g=4.5, m_sigma=1.0, n_fermions=1, N=400)
    assert strong["M"][0] < weak["M"][0]      # deeper well
    assert strong["bound"][0] < weak["bound"][0]   # lower ground-state energy


def test_pauli_filling_enlarges_and_shallows_binding():
    # filling more levels (Pauli) makes the soliton larger and shallows the binding per
    # particle -- the same direction the Thomas-Fermi Fermi-ball (fermi_ball.py) predicts
    tr = sc.degeneracy_trend(n_list=(1, 2, 4), g=4.0, g_v=0.0, N=500, R=10.0)
    assert all(tr["converged"])
    assert np.all(np.diff(tr["r_mean"]) > 0)               # grows with filling
    assert np.all(np.diff(tr["binding_per_particle"]) < 0)  # binding per particle shallows


def test_vector_channel_flattens_the_interior():
    # the torsion repulsion (vector channel, g_v>0) is what flattens the core into a box-like
    # drop: at every filling it raises M(0)/m0 (less deep) and lowers the central density
    # vs the scalar-only solver -- reproducing the fermi_ball direction the scalar model missed
    scalar = sc.degeneracy_trend(n_list=(1, 2, 4), g=4.0, g_v=0.0, N=500, R=10.0)
    vector = sc.degeneracy_trend(n_list=(1, 2, 4), g=4.0, g_v=2.5, N=500, R=10.0)
    assert all(vector["converged"])
    assert np.all(vector["interior_flatness"] > scalar["interior_flatness"])  # flatter core
    assert np.all(vector["rho0"] < scalar["rho0"])                            # less dense core


def test_crosscheck_direction_matches_fermi_ball():
    # the EOS shortcut and the full solver must agree on the DIRECTION: degeneracy makes the
    # drop less dense and more weakly bound
    from cartasis_sims import fermi_ball as fb
    shift = fb.degeneracy_shift()
    assert shift["n0_ratio"] < 1.0          # Pauli lowers the saturation density
    assert shift["binding_ratio"] < 1.0     # and shallows the binding
    tr = sc.degeneracy_trend(n_list=(1, 4), g=4.0, g_v=2.5, N=500, R=10.0)
    assert tr["binding_per_particle"][-1] < tr["binding_per_particle"][0]  # same direction
