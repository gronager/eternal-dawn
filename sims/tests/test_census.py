"""Tests for the supraverse census and our footprint."""

import math

from cartasis_sims import census as cs


def test_generation_pmf_is_geometric_and_bhu1_dominated():
    n, p = cs.generation_pmf(0.1)
    assert math.isclose(p[0], 0.9, rel_tol=1e-6)       # P(BHU1) = 1 - eps
    assert p[0] > p[1] > p[2]                           # decreasing tower


def test_chirality_is_5050():
    c = cs.chirality_split()
    assert math.isclose(c["matter"], 0.5) and math.isclose(c["antimatter"], 0.5)


def test_census_normalizes():
    rows = cs.census(n_max=8)
    assert math.isclose(sum(r[3] for r in rows), 1.0, rel_tol=1e-6)


def test_our_footprint_is_the_dominant_cell():
    fp = cs.our_footprint(0.1, 0.9)
    assert fp.generation == 1 and fp.chirality == "matter" and fp.family == "clean"
    # ~0.9 * 0.5 * 0.9 = 0.405 of all viable observers
    assert math.isclose(fp.fraction, 0.405, rel_tol=1e-6)
    assert cs.footprint_is_typical(0.1, 0.9)            # > 10% -> typical (Copernican)


def test_spin_does_not_change_generation():
    # the correction: a shared-axis null does NOT promote BHU1 -> BHU2
    assert not cs.spin_changes_generation()


def test_cpt_twin_has_equal_share():
    rows = cs.census(n_max=2)
    matter = sum(r[3] for r in rows if r[1] == "matter")
    anti = sum(r[3] for r in rows if r[1] == "antimatter")
    assert math.isclose(matter, anti, rel_tol=1e-9)     # identical by CPT
