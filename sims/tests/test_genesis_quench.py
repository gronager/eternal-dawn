"""The 2D baby-Skyrme genesis quench: topology, stability, condensation, Kibble-Zurek."""
import numpy as np
from cartasis_sims import genesis_quench as gq

STABLE = dict(kappa=5.0, mu2=0.01, dt=0.03, gamma=0.5)        # resolved baby-Skyrmion (size ~5)


def test_topological_charge_counts_winding():
    assert abs(abs(gq.total_charge(gq.baby_skyrmion(64, m=1))) - 1.0) < 0.2   # |Q|=1
    assert abs(abs(gq.total_charge(gq.baby_skyrmion(64, m=2))) - 2.0) < 0.3   # |Q|=2


def test_skyrmion_stable_not_collapsing():
    # a resolved Skyrmion keeps its winding under relaxation (does not fall through the lattice)
    sk = gq.baby_skyrmion(64, m=1, width=6.0)
    res = gq.run_quench(L=64, steps=1500, T0=0.0, T1=0.0, n_init=sk, hot_start=False, seed=1, **STABLE)
    assert abs(res["Q"][-1]) > 0.8                            # charge conserved, not unwound to 0


def test_quench_condenses_torsitons():
    # cooling a hot fluid drops the energy hugely and leaves a frozen relic of torsitons
    res = gq.run_quench(L=80, steps=3000, T0=0.9, T1=0.0, seed=2, record=6, **STABLE)
    assert res["energy"][-1] < 0.1 * res["energy"][0]        # condensed (energy collapses)
    assert res["content"][-1] > 0.5                          # relic torsitons survive
    assert abs(res["Q"][-1] - round(res["Q"][-1])) < 0.25    # net charge near-integer (topological)


def test_kibble_zurek_faster_freezes_more():
    rates, content = gq.quench_rate_scan([0.5, 2.0, 8.0], L=72, base_steps=6000, T0=0.9, seed=4, **STABLE)
    assert content[-1] > content[0]                          # faster quench -> more torsitons
