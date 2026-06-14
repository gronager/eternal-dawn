r"""A4/Z3 flavour: the Koide masses and tri-bimaximal mixing from ONE discrete symmetry.

A4 (tetrahedral, order 12) assigns the three generations to its triplet 3. Breaking A4 to different
residual subgroups in the two sectors produces masses and mixing together:

  * CHARGED LEPTONS (residual Z3): the mass matrix is Z3-CIRCULANT,
        M_e = circ(c0, c1, c1*) ,   eigenvalues  sqrt(m_k) = c0 + 2|c1| cos(phi + 2 pi k/3),
    which is EXACTLY the Koide/Foot form. Koide Q = 2/3 holds precisely at |c1|/c0 = 1/sqrt2 (the sqrt2
    amplitude = the circulant off/diagonal ratio). The phase phi = the Koide phase delta sets the
    hierarchy; every circulant is diagonalised by the Z3 Fourier ("magic") matrix.
  * NEUTRINOS (residual Z2 = mu-tau): the Majorana matrix is magic + mu-tau symmetric, FORM-diagonalised
    by the tri-bimaximal matrix -> the large PMNS angles (theta12: sin^2=1/3, theta23=45, theta13=0).

So one A4 gives BOTH: the Koide-form charged-lepton masses (Z3) and TBM mixing (Z3 vs Z2 misalignment).
Honest fine-print (the residual special points A4 does NOT force): the sqrt2 (exact 2/3) is a specific
singlet/doublet coupling ratio, and exact TBM has theta13=0 whereas the measured theta13~8.6deg is a
real correction (TBM broken at ~0.15, as `mixing.py` shows). A4 supplies the STRUCTURE that unifies
masses and mixing; the magnitude phase delta and the small TBM-breaking remain inputs.
"""
from __future__ import annotations

import numpy as np

W = np.exp(2j * np.pi / 3)
LEPTONS = np.array([0.51099895, 105.6583755, 1776.86])   # MeV


def circulant_sqrt_mass(c0, r, delta):
    """Z3-circulant sqrt-mass matrix circ(c0, c1, c1*) with c1 = c0 r e^{i delta}. Hermitian, real
    eigenvalues sqrt(m_k) = c0 + 2|c1| cos(delta + 2 pi k/3) -- the Koide form. r = |c1|/c0."""
    c1 = c0 * r * np.exp(1j * delta)
    return np.array([[c0, c1, np.conj(c1)],
                     [np.conj(c1), c0, c1],
                     [c1, np.conj(c1), c0]])


def magic_matrix():
    """The Z3 Fourier matrix that diagonalises ANY circulant (the charged-lepton diagonaliser)."""
    return np.array([[1, 1, 1], [1, W, W**2], [1, W**2, W]]) / np.sqrt(3)


def koide_Q(masses):
    m = np.abs(np.asarray(masses, float))
    s = np.sqrt(m)
    return float(m.sum() / s.sum() ** 2)


def charged_lepton_masses(c0, r, delta):
    """Squared eigenvalues of the circulant = the charged-lepton masses (up to overall scale c0^2)."""
    ev = np.linalg.eigvalsh(circulant_sqrt_mass(c0, r, delta))
    return np.sort(ev**2)


def fit_delta_to_leptons():
    """Find (scale, delta) so the circulant (r=1/sqrt2) reproduces the lepton masses; returns both."""
    from itertools import product
    sq = np.sqrt(LEPTONS)
    m0 = sq.mean()
    best = None
    for d in np.linspace(0, 2 * np.pi, 6000):
        pred = np.sort(np.linalg.eigvalsh(circulant_sqrt_mass(m0, 1 / np.sqrt(2), d)) ** 2)
        err = np.sum((np.log(pred) - np.log(np.sort(LEPTONS))) ** 2)
        if best is None or err < best[0]:
            best = (err, d, pred)
    return m0, best[1], best[2]


def neutrino_mutau_matrix(a, b, c):
    """A magic + mu-tau symmetric Majorana matrix: a*I + b*democratic + c*(mu-tau). TBM-diagonalised."""
    dem = np.ones((3, 3))
    mt = np.array([[0, 0, 0], [0, 1, -1], [0, -1, 1]], float)
    return a * np.eye(3) + b * dem + c * mt


def tbm():
    """The tri-bimaximal mixing matrix."""
    return np.array([[np.sqrt(2 / 3), np.sqrt(1 / 3), 0],
                     [-np.sqrt(1 / 6), np.sqrt(1 / 3), -np.sqrt(1 / 2)],
                     [-np.sqrt(1 / 6), np.sqrt(1 / 3), np.sqrt(1 / 2)]])


def report():
    print("A4/Z3: Koide masses (Z3) and TBM mixing (Z2) from one discrete symmetry\n")
    print("  CHARGED LEPTONS -- Z3-circulant sqrt-mass; Koide Q vs the off/diagonal ratio |c1|/c0:")
    for r in (0.4, 1 / np.sqrt(2), 0.9):
        Q = koide_Q(charged_lepton_masses(1.0, r, 0.0))   # delta=0: all eigenvalues positive
        tag = "  <- exact 2/3 (the sqrt2 amplitude)" if abs(r - 1 / np.sqrt(2)) < 1e-3 else ""
        print(f"    |c1|/c0 = {r:.3f}:  Koide Q = {Q:.4f}{tag}")
    m0, d, pred = fit_delta_to_leptons()
    print(f"\n  at r=1/sqrt2, the phase delta sets the hierarchy: delta={np.degrees(d):.0f} deg ->")
    print(f"    predicted lepton masses (MeV): {np.round(pred, 3)}")
    print(f"    observed                (MeV): {np.round(np.sort(LEPTONS), 3)}   (Koide Q={koide_Q(LEPTONS):.4f})")

    print("\n  NEUTRINOS -- magic + mu-tau Majorana -> tri-bimaximal:")
    Mn = neutrino_mutau_matrix(0.6, 0.25, 0.4)
    _, Un = np.linalg.eigh(Mn)
    print(f"    |U_nu - TBM| max = {np.abs(np.abs(Un) - np.abs(tbm())).max():.3f}  (form-diagonalised)")
    print("\n  => ONE A4: Z3 (charged) gives the Koide-form masses, Z2 (neutrino) gives TBM. The sqrt2")
    print("  (exact 2/3) and the small theta13 TBM-breaking are the residual inputs A4 does not force.")


if __name__ == "__main__":
    report()
