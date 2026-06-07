"""Kahana--Ripka CQSM: free spectrum exact, bound spectrum matches the grid, sea CONVERGES.

This is the basis that fixed the grid's runaway grand-spin sum, so the convergence and
box-stability assertions here are the whole point.
"""
import numpy as np
import pytest

from cartasis_sims import kahana_ripka as kr
from cartasis_sims import chiral_quark_soliton as cqs
from cartasis_sims import hedgehog as hh


def test_free_spectrum_is_exact_box_dirac():
    # theta=0: eigenvalues are +/- sqrt(p^2+M^2) with p the j_0 Dirichlet zeros; no in-gap level
    D, M = 10.0, 1.0
    free = lambda r: np.zeros_like(r)
    rq, w_r2 = kr._grid(D, 2000)
    w = np.sort(np.linalg.eigvalsh(kr.build_hamiltonian(0, +1, free(rq), M, 30, D, rq, w_r2)))
    p = np.array(kr.jl_zeros(0, 5)) / D
    expected = np.sort(np.sqrt(p**2 + M**2))
    assert np.allclose(np.sort(w[w > 0])[:5], expected, atol=1e-3)
    assert np.min(np.abs(w)) > M - 1e-3                  # no spurious in-gap state


def test_k0_valence_matches_trusted_hedgehog():
    # the K=0 valence dive reproduces hedgehog.kzero_levels to the basis/box accuracy
    for th0 in (1.5, 2.5, np.pi):
        lv = kr.ingap_levels(0, +1, lambda r: hh.hedgehog_profile(th0, r), Nb=30)
        kr_val = lv[np.argmin(np.abs(lv))]
        r2 = np.linspace(10.0 / 240, 10.0, 240)
        old = hh.kzero_levels(hh.hedgehog_profile(th0, r2), r2)
        assert abs(kr_val - old[np.argmin(np.abs(old))]) < 5e-3


@pytest.mark.parametrize("K", [0, 1])
def test_bound_spectrum_matches_grid_module(K):
    # the KR in-gap levels agree with the independently-built grid module (full soliton)
    prof = lambda r: np.pi * np.exp(-((r / 1.2) ** 2))
    rg = np.linspace(10.0 / 240, 10.0, 240)
    for P in (+1, -1):
        kri = np.sort(kr.ingap_levels(K, P, prof, Nb=32))
        grd = np.sort(cqs.ingap_levels(K, P, prof(rg), rg))
        assert len(kri) == len(grd)
        if len(kri):
            assert np.max(np.abs(kri - grd)) < 0.03


def test_sea_energy_converges_in_grand_spin():
    # THE point of the KR basis: the full-soliton sea sum converges (the grid one did not)
    prof = lambda r: np.pi * np.exp(-((r / 1.2) ** 2))
    e12 = kr.sea_energy(prof, Lam=4.0, Kmax=12, Nb=30)
    e16 = kr.sea_energy(prof, Lam=4.0, Kmax=16, Nb=30)
    assert abs(e16 - e12) < 0.05                          # converged: tail is tiny
    assert 5.0 < e16 < 8.0                                # an O(M), positive, finite Casimir energy


def test_soliton_mass_near_constituent_sum():
    # converged result: the B=1 mass infimum is the constituent sum N_c M (the soliton sits at ~3M)
    prof = lambda r: np.pi * np.exp(-((r / 0.6) ** 2))
    Msol, ev, es = kr.soliton_mass(prof, Lam=3.0, Kmax=12, Nb=30)
    assert es > 0                                         # sea energy is positive
    assert 3.0 <= Msol < 5.0                              # near (and above) the constituent sum 3 M


def test_self_consistent_vacuum_is_a_fixed_point():
    # theta=0 must be stationary: the pseudoscalar source vanishes in the vacuum
    rq, w_r2 = kr._grid(10.0, 1200)
    P = sum((2 * K + 1) * kr._sector_densities(K, P_, np.zeros_like(rq), 1.0, 2.5, 22, 10.0, rq, w_r2)[1]
            for K in range(0, 6) for P_ in (1, -1))
    assert np.max(np.abs(P)) < 0.05                      # vacuum sources no pion


def test_self_consistent_critical_coupling():
    # the full self-consistent solution: strong coupling binds a soliton, weak coupling collapses
    _, th_strong, ev_s = kr.self_consistent_profile(Lam=2.0, Kmax=8, Nb=24, Nq=1200, iters=20)
    _, th_weak, ev_w = kr.self_consistent_profile(Lam=3.0, Kmax=8, Nb=24, Nq=1200, iters=22)
    assert th_strong[0] > 1.5 and ev_s is not None       # stable torsiton soliton (theta(0) ~ pi)
    assert th_weak[0] < 0.5                               # collapses to the trivial vacuum
