"""Tests for the galaxy-spin de-confounding pipeline."""

import numpy as np

from cartasis_sims import galaxy_spins as gs


def test_equatorial_to_galactic_landmarks():
    # Galactic centre ~ (RA, Dec) = (266.405, -28.936) -> (l, b) ~ (0, 0)
    l, b = gs.equatorial_to_galactic(266.405, -28.936)
    assert abs(b) < 1.0
    assert min(l, abs(l - 360)) < 1.0
    # North Galactic pole ~ (192.859, 27.128) -> b ~ +90
    _, b2 = gs.equatorial_to_galactic(192.85948, 27.12825)
    assert b2 > 89.0


def _injected_catalogue(axis_lb, amplitude, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    g = rng.standard_normal((n, 3))
    vecs = g / np.linalg.norm(g, axis=1, keepdims=True)
    d = gs.lb_to_vec(*axis_lb)
    mu = vecs @ d
    p_cw = 0.5 * (1.0 + amplitude * mu)
    spin = np.where(rng.random(n) < p_cw, 1.0, -1.0)
    return vecs, spin, d


def test_fit_dipole_recovers_injected_axis():
    vecs, spin, d = _injected_catalogue((123.0, 45.0), amplitude=0.3)
    fit = gs.fit_dipole(vecs, spin)
    ang = np.degrees(np.arccos(abs(np.dot(fit["axis"], d))))
    assert ang < 8.0
    # 3|D| recovers the injected amplitude to ~20%
    assert abs(fit["amplitude_3D"] - 0.3) < 0.06


def test_null_pvalue_calibration():
    # random (isotropic) spins -> no significant dipole
    rng = np.random.default_rng(2)
    g = rng.standard_normal((8000, 3))
    vecs = g / np.linalg.norm(g, axis=1, keepdims=True)
    spin = rng.choice([-1.0, 1.0], size=8000)
    p = gs.dipole_null_pvalue(vecs, spin, n_perm=500, seed=3)
    assert p > 0.05
    # injected dipole -> significant
    vecs2, spin2, _ = _injected_catalogue((10.0, 20.0), amplitude=0.25)
    p2 = gs.dipole_null_pvalue(vecs2, spin2, n_perm=500, seed=4)
    assert p2 < 0.05
