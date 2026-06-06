"""The NJL-derived torsiton couplings: the Fierz ratio and the QCD-like check."""
import numpy as np

from cartasis_sims import njl_couplings as njl


def test_fierz_ratio_is_two():
    # the whole point: m_omega, g_v, m_sigma, lambda recombine to the Fierz channel ratio G_V/G_S=2
    c = njl.njl_couplings(g=4.0)
    assert np.isclose(c["GV_over_GS"], 2.0)
    assert np.isclose(c["g_v"], np.sqrt(3.0) * c["g"])
    assert np.isclose(c["m_omega"], np.sqrt(6.0) * c["M"])
    assert np.isclose(c["m_sigma"], 2.0 * c["M"])
    assert np.isclose(c["lam"], 2.0 * c["g"] ** 2)


def test_vector_heavier_than_scalar():
    # the resolution of the unbinding: the repulsive channel is shorter-range (heavier) than the
    # attractive one -- not the equal-mass naive assumption
    c = njl.njl_couplings(g=4.0)
    assert c["m_omega"] > c["m_sigma"]


def test_loop_gives_qcd_like_g():
    # a QCD-like cutoff (Lambda ~ 3-4 M) yields g = M/f_pi ~ 3-4 (M~330, f_pi~93 -> ~3.5)
    for LoM in (3.0, 4.0):
        c = njl.njl_couplings(Lambda_over_M=LoM)
        assert 2.0 < c["g"] < 5.0
        assert np.isclose(c["GV_over_GS"], 2.0)       # consistency holds at any cutoff
