"""A4/Z3: circulant -> Koide 2/3 at |c1|/c0=1/sqrt2; magic+mu-tau -> tribimaximal."""
import numpy as np

from cartasis_sims import a4_flavor as a4


def test_circulant_gives_koide_two_thirds_at_sqrt2():
    # Koide Q = 2/3 exactly at the circulant off/diagonal ratio 1/sqrt2 (delta=0 keeps eigenvalues +)
    Q = a4.koide_Q(a4.charged_lepton_masses(1.0, 1 / np.sqrt(2), 0.0))
    assert abs(Q - 2 / 3) < 1e-3
    # away from 1/sqrt2 it departs from 2/3
    assert abs(a4.koide_Q(a4.charged_lepton_masses(1.0, 0.4, 0.0)) - 2 / 3) > 0.05


def test_circulant_phase_reproduces_lepton_masses():
    m0, d, pred = a4.fit_delta_to_leptons()
    assert np.allclose(np.sort(pred), np.sort(a4.LEPTONS), rtol=0.02)


def test_magic_diagonalises_circulant():
    F = a4.magic_matrix()
    M = a4.circulant_sqrt_mass(1.0, 1 / np.sqrt(2), 0.3)
    D = F.conj().T @ M @ F
    assert np.abs(D - np.diag(np.diag(D))).max() < 1e-9      # off-diagonal vanishes


def test_mutau_neutrino_gives_tribimaximal():
    _, Un = np.linalg.eigh(a4.neutrino_mutau_matrix(0.6, 0.25, 0.4))
    assert np.abs(np.abs(Un) - np.abs(a4.tbm())).max() < 1e-6
