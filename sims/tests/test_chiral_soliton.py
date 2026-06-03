"""Tests for the chiral (symmetry-breaking) self-consistent soliton."""

import math

import numpy as np

from cartasis_sims import chiral_soliton as ch


def test_chiral_loop_converges_and_binds():
    out = ch.solve_chiral(g=4.0, lam=4.0, n_fermions=2, eps_break=0.0, N=400)
    assert out["converged"]
    assert len(out["bound"]) >= 1
    assert out["f_pi"] == 1.0


def test_mass_is_configurational():
    # the observable soliton mass is overwhelmingly FIELD + BAG energy, not constituent
    out = ch.solve_chiral(g=5.0, lam=4.0, n_fermions=4, eps_break=0.0, N=500)
    field_frac = (out["grad_energy"] + out["pot_energy"]) / out["mass"]
    assert field_frac > 0.8                       # like the proton (~99% binding)
    assert out["mass"] > 10 * out["constituent_mass"]   # mass >> constituent gv


def test_mass_scales_with_condensate():
    # mass is created by the condensate: mass ~ linear in v, -> 0 as v -> 0
    m1 = ch.solve_chiral(v=1.0, g=5.0, lam=4.0, n_fermions=4, eps_break=0.0, N=400)["mass"]
    m2 = ch.solve_chiral(v=2.0, g=5.0, lam=4.0, n_fermions=4, eps_break=0.0,
                         R=6.0, N=400)["mass"]
    assert math.isclose(m2 / m1, 2.0, rel_tol=0.05)   # mass proportional to v


def test_s_estimate_is_flagged_proxy():
    # the chiral model gives a real f_pi, but the M_V proxy is the fermion gap (wrong
    # scale), so the S estimate is NOT trustworthy -- guard that it is not mistaken
    # for a delivered number by checking it sits in the unphysical regime here.
    out = ch.solve_chiral(g=5.0, lam=4.0, n_fermions=4, eps_break=0.0, N=400)
    s = ch.s_estimate(out)
    assert s["f_pi"] == 1.0
    assert s["S_proxy"] > 1.0                      # unphysically large -> a proxy, not S
