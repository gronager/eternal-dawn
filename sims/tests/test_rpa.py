"""Tests for the RPA-dressed S comparison (torsion vs QCD-like couplings)."""

import math

from cartasis_sims import rpa


def test_leading_limit_is_Nc_over_6pi():
    # at zero coupling the RPA reduces to the leading constituent-loop value
    S = rpa.s_rpa(0.0, 0.0)
    assert math.isclose(S, 3 / (6 * math.pi), rel_tol=1e-6)


def test_walking_loops_equalise():
    # large Lambda/M (small chiral breaking) -> Pi_A/Pi_V -> 1
    near = rpa.compare(3.0)["PiA_over_PiV"]
    deep = rpa.compare(1000.0)["PiA_over_PiV"]
    assert deep > near
    assert deep > 0.9


def test_equal_couplings_beat_qcd_in_walking_limit():
    # the torsion (G_A=G_V) advantage grows as the sector walks
    near = rpa.compare(3.0)["ratio"]
    deep = rpa.compare(1000.0)["ratio"]
    assert deep < near                       # torsion gets relatively better
    assert deep < 0.5                        # ~3x smaller S than QCD-like deep in walking


def test_torsion_S_below_qcd_S():
    c = rpa.compare(100.0)
    assert c["S_torsion"] < c["S_qcd"]
