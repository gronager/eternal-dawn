"""The hedgehog: the grand-spin K=0 valence level dives into the gap and binds."""
import numpy as np

from cartasis_sims import hedgehog as hh


def test_free_limit_has_no_bound_level():
    # theta=0 is the free quark: no level inside the mass gap
    r = np.linspace(10.0 / 240, 10.0, 240)
    assert len(hh.kzero_levels(hh.hedgehog_profile(0.0, r), r)) == 0


def test_kzero_level_dives_and_binds_deeply():
    # as the chiral angle grows, a K=0 level dives from the continuum deep into the gap
    r = np.linspace(10.0 / 240, 10.0, 240)
    deep = hh.valence_energy(np.pi, r)
    assert not np.isnan(deep) and abs(deep) < 0.5          # full hedgehog: deeply bound (~0.2 M)


def test_diving_is_monotonic_in_chiral_amplitude():
    # the valence level dives deeper as theta0 increases (the soliton grows)
    r = np.linspace(10.0 / 240, 10.0, 240)
    es = [hh.valence_energy(t, r) for t in (1.5, 2.0, 2.5, 3.0)]
    assert all(not np.isnan(e) for e in es)
    assert es[0] > es[1] > es[2] > es[3]                   # monotonic dive
