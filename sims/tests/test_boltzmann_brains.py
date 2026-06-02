"""Tests for the Boltzmann-brain / arrow-priority resolution (three towers)."""

from cartasis_sims import boltzmann_brains as bb


def test_three_tower_inner_exponents():
    t = bb.three_towers()
    # brain ~69, seed ~84, bare Big Bang 123 -- the chapter's three numbers
    assert 68.0 < t["brain"] < 71.0
    assert 83.0 < t["seed"] < 85.0
    assert t["big_bang"] == 123.0


def test_brain_is_cheapest_seed_in_the_middle():
    t = bb.three_towers()
    assert t["brain"] < t["seed"] < t["big_bang"]


def test_mass_energy_dominates_configuration_per_particle():
    # per nucleon, the de Sitter mass toll (~5e42) crushes the ~1-bit config toll,
    # so a brain can never be made expensive by its organization
    assert bb.per_particle_mass_vs_config() > 1e40


def test_single_tower_rescues_cannot_close_the_gap():
    assert bb.single_tower_rescues_fail()


def test_finite_dark_energy_starves_interior_brains():
    # brain recurrence ~10^(10^69) yr >> dark-energy lifetime ~1e136 yr
    assert bb.brain_recurrence_time_years() > 1e69
    assert bb.finite_dark_energy_beats_brains()
