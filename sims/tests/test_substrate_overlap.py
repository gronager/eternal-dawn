"""Tests for the substrate-overlap pass (the generation hierarchy from the sharp substrate)."""

import numpy as np

from cartasis_sims import substrate_overlap as su


def test_healing_length_is_derived_from_well():
    # xi = width / sqrt(depth) -- a derived scale, not a fitted one
    assert np.isclose(su.healing_length(depth=4.0, width=1.0), 0.5)


def test_substrate_spread_is_orders_of_magnitude():
    # the substrate overlap gives a steep geometric hierarchy (~10^3), not ~3
    m = su.substrate_overlap_masses(depth=6.0)
    assert su.spread(m) > 300


def test_substrate_beats_broad_self_condensate():
    c = su.compare_to_self_condensate(depth=6.0)
    # the sharp substrate sharpens the spread far beyond the broad self-condensate
    assert c["substrate_spread"] > 100 * c["broad_spread"]
    assert c["broad_spread"] < 5


def test_masses_are_ascending_and_normalised():
    m = su.substrate_overlap_masses(depth=6.0)
    assert np.all(np.diff(m) > 0)        # gen1 < gen2 < gen3
    assert np.isclose(m.max(), 1.0)


def test_substrate_reaches_the_10e3_ballpark():
    # the substrate overlap reaches a ~10^3 spread at a natural O(1-10) depth -- the right
    # order of magnitude for the generation hierarchy (matches down-type ~889, within ~3x of
    # leptons ~3477). It does NOT reach up-type's ~80000 with this derived width: that extra
    # steepness is the owed inter-tower/dynamics refinement, recorded honestly here.
    d, s = su.depth_for_spread(su.OBSERVED_SPREAD["down-quark"])
    assert 3.0 <= d <= 30.0
    assert 0.5 < s / su.OBSERVED_SPREAD["down-quark"] < 2.0     # down-type matched
    best = max(su.spread(su.substrate_overlap_masses(dd)) for dd in (4.0, 6.0, 8.0))
    assert best > 300          # reaches the 10^3 ballpark (vs the broad overlap's ~3)
