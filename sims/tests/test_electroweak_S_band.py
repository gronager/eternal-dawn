"""L3 electroweak S: the honest band, and the relativistic cutoff tell."""
import numpy as np

from cartasis_sims import electroweak_S as ew


def test_leading_S_in_graveyard():
    assert np.isclose(ew.s_leading(Nc=3), 3 / (6 * np.pi))
    assert ew.s_leading(Nc=3) > ew.S_LEP_BOUND          # ~0.16 > 0.1: the graveyard


def test_walking_band_spans_graveyard_to_escape():
    band = ew.s_band(MV_over_fpi=(8.0, 16.0))
    # QCD-like (M_V/f_pi ~ 8) sits in the graveyard; deep walking (~16) escapes below the LEP bound
    assert band["S_max"] > ew.S_LEP_BOUND
    assert band["S_min"] < ew.S_LEP_BOUND
    # the escape threshold is well above QCD's 8.3 -- walking is required, not optional
    assert band["MV_over_fpi_for_S0p1"] > band["qcd_MV_over_fpi"]
    assert 12.0 < band["MV_over_fpi_for_S0p1"] < 15.0


def test_walking_form_matches_known_points():
    # QCD-like M_V/f_pi = 8.3 -> S ~ 0.25-0.3 (graveyard); M_V/f_pi = 15 -> S < 0.1 (escape)
    assert 0.2 < ew.s_walking(1 / 8.3) < 0.35
    assert ew.s_walking(1 / 15.0) < ew.S_LEP_BOUND


def test_first_weinberg_sum_rule_diverges_on_the_loop():
    # the relativistic tell: int(rho_V - rho_A) ds on the bare continuum grows ~log(cutoff) --
    # it does NOT converge, so a cutoff-free escape cannot come from the loop alone
    cutoffs = np.array([5.0, 50.0, 500.0])           # in units of M
    running = ew.wsr_continuum_divergence(1.0, cutoffs)
    assert running[1] > running[0] and running[2] > running[1]   # keeps growing with the cutoff
    # log growth: roughly constant increments per decade, not saturating
    assert (running[2] - running[1]) > 0.3 * (running[1] - running[0])
