"""A4 invariance of the Hehl-Datta flavour structure: group relations, Schur, residual -> circulant/mu-tau."""
import numpy as np

from cartasis_sims import a4_symmetry as a4s


def test_generators_close_to_A4():
    rel = a4s.group_relations()
    assert rel["S^2=I"] and rel["T^3=I"] and rel["(ST)^3=I"]
    assert rel["order"] == 12                      # <S,T> = A4


def test_schur_forces_degenerate_mass():
    # commutant of the irreducible triplet is 1-dimensional (scalars) -> A4-invariant M = m I
    nd, _ = a4s.schur_mass()
    assert nd == 1


def test_four_fermion_singlet_invariant():
    assert a4s.four_fermion_invariant()


def test_residual_symmetries_fix_matrix_forms():
    cd, iscirc, _ = a4s.circulant_from_Z3()
    assert cd == 3 and iscirc                      # commutant of Z3 = circulants (Koide)
    assert a4s.mutau_from_Z2() == 5                # commutant of Z2 = mu-tau forms (TBM)


def test_generators_unitary():
    for g in (a4s.S(), a4s.T(), a4s.P()):
        assert np.allclose(g.conj().T @ g, np.eye(3))
