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


# --- OGU mass distribution -------------------------------------------------
def test_mass_density_is_declining_power_law():
    M = np.logspace(0.5, 2.5, 200)        # above the birth scale
    n = pop.ogu_mass_density(M, g=2.0, M0=1.0, M_max=1e4)
    assert np.all(np.diff(n) < 0)         # many small, few large
    # log-log slope near -g away from the cutoffs
    lo, hi = 40, 120
    slope = np.polyfit(np.log(M[lo:hi]), np.log(n[lo:hi]), 1)[0]
    assert -2.3 < slope < -1.7


def test_observer_mass_peaks_at_viability_edge():
    M = np.logspace(0, 3, 400)
    no = pop.observer_mass_density(M, g=2.0, M_vis=10.0)
    # declining power law gated below M_vis -> peak sits right at the edge
    assert abs(M[np.argmax(no)] - 10.0) / 10.0 < 0.1
    assert np.all(no[M < 10.0] == 0.0)


# --- Spin / viability ------------------------------------------------------
def test_spin_birth_peaks_at_zero():
    w = np.linspace(-4, 4, 401)
    p = pop.spin_birth_pdf(w, sigma=1.0)
    assert abs(w[np.argmax(p)]) < 1e-6    # low spin likeliest


def test_observer_spin_peaks_above_threshold():
    w = np.linspace(0, 5, 2000)
    wmin = 0.6
    p = pop.spin_observer_pdf(w, sigma=1.0, omega_min=wmin, prod_power=1.0)
    wpeak = w[np.argmax(p)]
    assert wpeak > wmin                   # peak is past the sterile band
    assert wpeak < 2.0                    # but still low spin (near the edge)


def test_sterile_fraction_increases_with_threshold():
    assert pop.sterile_fraction(1.0, 1.0) > pop.sterile_fraction(1.0, 0.3)
    assert 0.0 < pop.sterile_fraction(1.0, 0.5) < 1.0


def test_purity_proportional_to_spin():
    assert np.isclose(pop.purity_from_spin(2.0, C=3.0, T=4.0),
                      pop.purity_from_spin(1.0, C=3.0, T=4.0) * 2.0)
