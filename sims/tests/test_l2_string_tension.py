"""L2: sigma = 2 pi v^2 -- the condensate sets the string tension (one substrate does both)."""
import numpy as np
import pytest

from cartasis_sims import l2_string_tension as l2
from cartasis_sims import chiral_soliton as cs


def test_bps_string_tension_is_2pi_v2():
    assert np.isclose(l2.bps_string_tension(1.0), 2 * np.pi)
    assert np.isclose(l2.bps_string_tension(246.0), 2 * np.pi * 246.0 ** 2)


def test_sigma_fpi_ratio_predicts_2pi():
    # the dimensionless lattice observable: sigma = 2 pi f_pi^2  ->  ratio = 2 pi
    assert np.isclose(l2.sigma_fpi_ratio(2 * np.pi * 0.09 ** 2, 0.09), 2 * np.pi)
    # a QCD-like vacuum (sigma ~ (440 MeV)^2, f_pi ~ 93 MeV) is NOT 2 pi -> the prediction is
    # non-trivial / falsifiable, not a generic identity
    assert l2.sigma_fpi_ratio(0.440 ** 2, 0.093) > 3 * (2 * np.pi)


def test_ano_vortex_saturates_bps_bound():
    # the testable known answer: a winding-1 vortex AT BPS (beta=1) has tension exactly 2 pi v^2
    assert np.isclose(l2.ano_vortex_tension(1.0), 1.0, atol=2e-3)


def test_ano_vortex_type_i_below_type_ii_above_bps():
    eps = {b: l2.ano_vortex_tension(b) for b in (0.5, 1.0, 2.0)}
    assert eps[0.5] < eps[1.0] < eps[2.0]        # monotonic in beta
    assert eps[0.5] < 1.0                          # type I (beta<1): below the bound
    assert eps[2.0] > 1.0                          # type II (beta>1): above the bound


def test_prediction_from_chiral_soliton():
    # pull v = f_pi from the chiral soliton and form the L2 prediction sigma = 2 pi v^2
    out = l2.predicted_sigma_from_soliton(cs.solve_chiral)
    assert np.isclose(out["v"], 1.0)               # default soliton has v = f_pi = 1
    assert np.isclose(out["sigma_bps"], 2 * np.pi)
