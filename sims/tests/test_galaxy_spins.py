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


def _data(fname):
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "..", "..", "data", "galaxy_spins", fname)


def test_load_both_real_catalogue_schemas():
    l1, b1, s1 = gs.load_spin_catalogue(_data("gz1_raw_debiased.csv"))
    l2, b2, s2 = gs.load_spin_catalogue(_data("iye2019.csv"))
    assert len(s1) > 30000 and set(np.unique(s1)) <= {-1.0, 1.0}      # GZ1 clean
    assert 450 < len(s2) < 600 and set(np.unique(s2)) <= {-1.0, 1.0}  # Iye 530
    assert np.all(np.abs(b1) <= 90.0) and np.all(np.abs(b2) <= 90.0)


def test_independent_methods_disagree_on_asymmetry():
    # GZ1 (crowd) shows a highly SIGNIFICANT CW excess (perception bias amplified by
    # its huge N); Iye (expert) is consistent with isotropy. The contrast is the
    # significance, not the raw fraction.
    from scipy.stats import binomtest
    _, _, s_gz1 = gs.load_spin_catalogue(_data("gz1_raw_debiased.csv"))
    _, _, s_iye = gs.load_spin_catalogue(_data("iye2019.csv"))
    p_gz1 = binomtest(int((s_gz1 > 0).sum()), len(s_gz1), 0.5).pvalue
    p_iye = binomtest(int((s_iye > 0).sum()), len(s_iye), 0.5).pvalue
    assert p_gz1 < 1e-6      # GZ1: a 'significant' asymmetry...
    assert p_iye > 0.05      # ...that the independent method does not reproduce


def test_both_axes_track_the_galactic_pole_not_cmb():
    # the honest systematic signature: best-fit axes sit near |b|~90, not the CMB axis
    for fname in ("gz1_raw_debiased.csv", "iye2019.csv"):
        l, b, spin = gs.load_spin_catalogue(_data(fname))
        vecs = np.array([gs.lb_to_vec(li, bi) for li, bi in zip(l, b)])
        fit = gs.fit_dipole(vecs, spin)
        assert abs(fit["b"]) > 70.0                # within ~20 deg of the Galactic pole
