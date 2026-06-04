"""Tests for the S-budget generation cap."""

import math

from cartasis_sims import generations_cap as gc


def test_leading_order_forbids_our_existence():
    # leading-order S (~0.16) gives a cap of zero -> walking is mandatory
    assert gc.n_max(gc.S_LEADING) == 0


def test_cap_three_window():
    lo, hi = gc.s_per_gen_for_cap(3)
    assert gc.n_max(0.5 * (lo + hi)) == 3
    assert lo < 0.033 <= hi


def test_cap_is_steep_in_s():
    # the 1/x sensitivity: halving S_per_gen roughly doubles the cap
    assert gc.n_max(0.05) == 2
    assert gc.n_max(0.025) == 4


def test_walking_is_required_for_three():
    w = gc.requires_walking()
    assert w["leading_cap"] == 0
    assert w["walk_factor_needed"] > 4         # need ~5x walk-down to allow 3
