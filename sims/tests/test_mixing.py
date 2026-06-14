"""Transitions/oscillations: 3 knots -> 2 gaps, 2+1 spectrum; PMNS mixed (TBM), CKM aligned."""
import numpy as np

from cartasis_sims import mixing as mx


def test_two_gaps_three_knots_2plus1():
    g = mx.mass_gap_structure()
    assert g["n_frequencies"] == 2 and g["n_knots"] == 3
    assert g["ratio"] > 20                       # 2+1: one gap >> the other (two close + one far)


def test_pmns_near_tribimaximal_ckm_aligned():
    pmns = mx.mixing_matrix(*mx.PMNS_DEG)
    ckm = mx.mixing_matrix(*mx.CKM_DEG)
    tbm = mx.tribimaximal()
    assert np.abs(np.abs(pmns) - np.abs(tbm)).max() < 0.2     # PMNS near tri-bimaximal
    assert mx.offdiagonal_power(pmns) > 10 * mx.offdiagonal_power(ckm)  # leptons mixed, quarks aligned


def test_ckm_nearly_identity():
    ckm = mx.mixing_matrix(*mx.CKM_DEG)
    assert np.abs(ckm - np.eye(3)).max() < 0.25  # quarks nearly aligned (mass ~ flavour basis)
