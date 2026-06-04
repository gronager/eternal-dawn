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
