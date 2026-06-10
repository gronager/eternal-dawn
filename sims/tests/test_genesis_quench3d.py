"""The 3D Skyrme genesis: baryon-number topology, Skyrmion stability, and quench condensation."""
import numpy as np
from cartasis_sims import genesis_quench3d as g3


def test_hedgehog_is_one_baryon():
    # a hedgehog has |B|=1 (recovered from below by finite-difference discretization)
    assert 0.7 < abs(g3.total_baryon(g3.hedgehog_skyrmion(32, width=4.0))) < 1.15


def test_skyrmion_stable_not_collapsing():
    # the Skyrme stabiliser holds a resolved 3D baryon together under relaxation (B conserved)
    sk = g3.hedgehog_skyrmion(32, width=4.5)
    res = g3.run_quench3d(L=32, steps=300, T0=0.0, T1=0.0, n_init=sk, hot_start=False, seed=1, record=3)
    assert abs(res["B"][-1]) > 0.7                       # not unwound to vacuum


def test_quench_condenses_baryons():
    # cooling a hot SU(2) chiral fluid drops the energy and leaves a relic of baryons (torsitons)
    res = g3.run_quench3d(L=28, steps=700, T0=0.85, T1=0.0, seed=2, record=4)
    assert res["energy"][-1] < 0.2 * res["energy"][0]    # condensed (energy collapses)
    assert res["content"][-1] > 0.5                      # baryons survive
    assert abs(res["B"][-1]) <= res["content"][-1] + 1e-9  # net |B| <= total winding (+/- pairs)
