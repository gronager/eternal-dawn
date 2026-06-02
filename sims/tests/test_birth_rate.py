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
    assert math.isclose(br.instanton_action_value(br.LAMBDA_CRIT), 1.43, abs_tol=0.05)


def test_separation_grows_as_lambda_shrinks():
    assert math.isclose(br.ogu_separation_hubble(1.0), 1.0, rel_tol=1e-9)
    assert br.ogu_separation_hubble(1e-20) > br.ogu_separation_hubble(0.24) > 1.0


def test_regime_labels_flip_at_threshold():
    assert "packed" in br.regime(0.5)
    assert "dilute" in br.regime(1e-10)


def test_minimal_bounce_seed_mass():
    # self-gravitating at rho_C (R_s = R): M ~ c^3/(2 G Omega_bounce) ~ 1e15 kg
    M = br.minimal_seed_mass(1.0e50)
    assert 1e14 < M < 1e16


def test_instanton_action_is_astronomically_above_threshold():
    # I = 2 pi M_seed c^2/(hbar H_Lambda) ~ 1e84 >> I_crit ~ 1.4  -> dilute, robustly
    I = br.seed_instanton_action()
    assert I > 1e80
    assert I / br.instanton_action_threshold() > 1e80


def test_boltzmann_action_linear_in_mass():
    assert math.isclose(br.de_sitter_boltzmann_action(2e14),
                        2.0 * br.de_sitter_boltzmann_action(1e14), rel_tol=1e-9)


def test_threshold_seed_is_sub_planckian_hence_cannot_bounce():
    # the seed that would give I_crit is far below the Planck mass -> no horizon,
    # cannot bounce, so the realised I is always >> I_crit
    assert br.seed_mass_for_threshold() < 1e-40       # << M_Planck ~ 2e-8 kg
