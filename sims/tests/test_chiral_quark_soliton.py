"""The CQSM grand-spin machinery: angular exactness, free spectrum, K=0 match, small-amplitude f_pi.

The validated foundation for the torsiton mass. The absolute sea Casimir energy of the full soliton
is NOT yet trustworthy on this grid (see module docstring) and is deliberately not asserted here.
"""
import numpy as np
import pytest

from cartasis_sims import chiral_quark_soliton as cqs
from cartasis_sims import hedgehog as hh


@pytest.mark.parametrize("K", [0, 1, 2])
def test_tau_rhat_squares_to_identity(K):
    # (tau.rhat)^2 = 1: the decisive correctness check on the grand-spin angular matrix elements
    S = cqs.grandspin_states(K)
    T = np.array([[cqs.tau_rhat(K, s2, s1) for s1 in S] for s2 in S])
    assert np.allclose(T @ T, np.eye(len(S)), atol=1e-9)


@pytest.mark.parametrize("K", [0, 1, 2])
def test_free_spectrum_has_no_in_gap_and_is_parity_degenerate(K):
    # theta=0 is the free quark: no level inside the mass gap, and the two parity sectors degenerate
    r = np.linspace(12.0 / 240, 12.0, 240)
    th = np.zeros_like(r)
    assert len(cqs.ingap_levels(K, +1, th, r)) == 0
    assert len(cqs.ingap_levels(K, -1, th, r)) == 0
    wp = np.sort(np.linalg.eigvalsh(cqs.build_hamiltonian(K, +1, th, r)))
    wm = np.sort(np.linalg.eigvalsh(cqs.build_hamiltonian(K, -1, th, r)))
    assert np.allclose(wp, wm, atol=1e-6)                 # free Dirac: parity-degenerate


def test_k0_sector_matches_trusted_hedgehog():
    # the K=0 (parity +) valence dive reproduces the independently-trusted hedgehog.kzero_levels
    r = np.linspace(10.0 / 240, 10.0, 240)
    for th0 in (1.5, 2.5, np.pi):
        th = hh.hedgehog_profile(th0, r)
        new = cqs.ingap_levels(0, +1, th, r)
        old = hh.kzero_levels(th, r)
        nv = new[np.argmin(np.abs(new))]
        ov = old[np.argmin(np.abs(old))]
        assert abs(nv - ov) < 2e-3


def test_small_amplitude_sea_is_a_clean_gradient_term():
    # for a REGULAR small-amplitude profile (theta(0)=0) E_sea ~ theta0^2 (a gradient term, no
    # spurious mass term) with an effective f_pi in the expected ~0.4 M band (proper-time scheme)
    Nr, Rb = 200, 10.0
    r = np.linspace(Rb / Nr, Rb, Nr)
    g = (r / 2.0) * np.exp(-((r / 2.0) ** 2))            # vanishes at origin AND infinity: regular
    dg = np.gradient(g, r)
    geom = np.trapezoid(4 * np.pi * r**2 * (dg**2 + 2 * g**2 / r**2), r)
    es = {th0: cqs.sea_energy(th0 * g, r, Lam=4.0, Kmax=12) for th0 in (0.15, 0.30)}
    # quadratic scaling: E_sea(0.30)/E_sea(0.15) ~ 4
    assert 3.5 < es[0.30] / es[0.15] < 4.5
    fpi2 = 2 * (es[0.30] / 0.30**2) / geom
    assert 0.15 < fpi2 < 0.30                            # f_pi ~ 0.39 - 0.55 M (proper-time band)
