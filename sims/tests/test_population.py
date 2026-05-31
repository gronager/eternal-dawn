"""Tests for the supraverse population / generation-depth model."""

import numpy as np

from cartasis_sims import population as pop


def test_pmf_normalized():
    p = pop.generation_pmf(m=1.5, D=8)
    assert np.isclose(p.sum(), 1.0)


def test_supercritical_pushes_deep():
    # m>1: typical universe sits near the truncation depth
    D = 10
    en = pop.expected_generation(m=2.0, D=D)
    assert en > D / 2


def test_subcritical_stays_shallow():
    en = pop.expected_generation(m=0.4, D=20)
    assert en < 1.5            # lineage dies out -> shallow


def test_structure_decay_pulls_observers_shallower():
    deep = pop.observer_weighted_pmf(m=2.0, D=12, structure_decay=1.0)
    shal = pop.observer_weighted_pmf(m=2.0, D=12, structure_decay=0.4)
    n = np.arange(13)
    assert (n * shal).sum() < (n * deep).sum()


def test_dm_de_rules_out_ogu():
    # the decisive anchor: we have DM and DE, so we are not an OGU
    assert pop.we_are_ogu(has_dark_matter=False, has_dark_energy=False) is True
    assert pop.we_are_ogu(has_dark_matter=True, has_dark_energy=True) is False
    assert pop.min_generation_from_observations(True, True) == 1


def test_posterior_zeros_ogu():
    p = pop.depth_posterior(m=2.0, D=10, n_min=1)
    assert p[0] == 0.0
    assert np.isclose(p.sum(), 1.0)


def test_survival_truncation_monotone():
    # higher per-generation viability -> deeper allowed tree
    assert pop.survival_truncation_depth(0.9) > pop.survival_truncation_depth(0.5)
