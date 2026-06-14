"""All 12 as 4 towers: leptons A4-exact (A=sqrt2), quarks/neutrinos broken."""
import numpy as np

from cartasis_sims import twelve_masses as tm


def test_leptons_are_a4_exact():
    c0, A, d, Q = tm.fit_circulant(tm.TOWERS["charged leptons"])
    assert abs(A - np.sqrt(2)) < 0.02          # amplitude = sqrt2 (Koide 2/3) exactly
    assert abs(Q - 2 / 3) < 1e-3


def test_quarks_and_neutrinos_are_broken():
    for name in ("up quarks", "down quarks", "neutrinos (NO)"):
        _, A, _, _ = tm.fit_circulant(tm.TOWERS[name])
        assert abs(A - np.sqrt(2)) > 0.05      # amplitude deviates from sqrt2 (A4 broken)


def test_all_towers_fit_circulant():
    # the circulant reproduces each tower's masses (3 params, exact structural fit)
    for name, m in tm.TOWERS.items():
        c0, A, d, Q = tm.fit_circulant(m)
        pred = np.sort([(c0 * (1 + A * np.cos(np.radians(d) + 2 * np.pi * k / 3))) ** 2 for k in range(3)])
        assert np.allclose(pred, np.sort(m), rtol=0.05)
