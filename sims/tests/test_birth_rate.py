"""Tests for pinning beta as a de Sitter nucleation rate."""

import math

from cartasis_sims import birth_rate as br


def test_hubble_lambda_from_dark_energy():
    # H_Lambda = H_0 sqrt(Omega_Lambda) ~ 1.8e-18 s^-1
    assert math.isclose(br.hubble_lambda(), 1.8e-18, rel_tol=0.05)


def test_prefactor_pinned_by_dark_energy():
    # H_Lambda^4/c^3 ~ 4e-97 m^-3 s^-1, fixed (not free)
    assert math.isclose(br.nucleation_prefactor(), 4.0e-97, rel_tol=0.15)


def test_beta_linear_in_lambda():
    assert math.isclose(br.beta(0.5), 0.5 * br.nucleation_prefactor(), rel_tol=1e-9)
    assert br.beta(1e-20) < br.beta(1e-3)


def test_percolation_threshold_and_action():
    assert br.percolates(0.5) and not br.percolates(0.01)
    assert math.isclose(br.instanton_action(br.LAMBDA_CRIT), 1.43, abs_tol=0.05)


def test_separation_grows_as_lambda_shrinks():
    assert math.isclose(br.ogu_separation_hubble(1.0), 1.0, rel_tol=1e-9)
    assert br.ogu_separation_hubble(1e-20) > br.ogu_separation_hubble(0.24) > 1.0


def test_regime_labels_flip_at_threshold():
    assert "packed" in br.regime(0.5)
    assert "dilute" in br.regime(1e-10)
