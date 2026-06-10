"""Confinement (Horizon 3a): quarks confined, baryons/mesons/leptons free -- the species distinction."""
import numpy as np
from cartasis_sims import genesis_confine as gc


def test_only_the_quark_is_confined():
    sigma = 0.08
    assert gc.is_confined("quark", sigma) is True          # no finite free-quark state
    for free in ("lepton", "meson", "baryon"):
        assert gc.is_confined(free, sigma) is False         # colour singlets / colourless are free
    # the quark energy grows linearly with the box (the string), the others are box-independent
    grow = gc.species_energy("quark", sigma, box=100) - gc.species_energy("quark", sigma, box=10)
    assert np.isclose(grow, sigma * 90)                    # energy grows as sigma x (box) -- the string
    for free in ("lepton", "meson", "baryon"):
        assert np.isclose(gc.species_energy(free, sigma, box=10), gc.species_energy(free, sigma, box=100))


def test_meson_has_a_finite_bound_size():
    assert np.isclose(gc.meson_size(0.08, core=0.5), np.sqrt(0.5 / 0.08))
    assert gc.string_potential(gc.meson_size(0.08), 0.08) < gc.string_potential(0.1, 0.08)  # it's the min


def test_binding_binds_quarks_not_leptons():
    res = gc.run_binding(n_quark_pairs=8, n_lepton=6, steps=2500, seed=2)
    assert res["bound_fraction"] > 0.7                      # coloured knots pulled into hadrons
