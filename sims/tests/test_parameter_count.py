"""Tests for the SM-vs-ED particle-sector parameter count."""

from cartasis_sims import parameter_count as pc


def test_sm_count_with_and_without_neutrinos():
    # ~19 without neutrino masses, ~26 with
    assert pc.sm_count(with_neutrinos=False) == 19
    assert pc.sm_count(with_neutrinos=True) == 26


def test_ed_adds_zero_new_fundamental_parameters_ideal():
    assert pc.ed_new_free_parameters_ideal() == 0


def test_base_constants_are_shared_not_new():
    # the base constants are G, hbar, c -- inherited from Part I / all of physics
    base = pc.ed_base_constants()
    assert base == ["G", "hbar", "c"]
    # the condensate scale is NOT among them (it is generated, not inserted)
    assert not any("f_pi" in b or "v" == b for b in base)


def test_replaced_count_matches_full_sm():
    # the programme aims to derive every SM parameter, neutrinos included
    assert pc.replaced_count() == pc.sm_count(with_neutrinos=True)


def test_residuals_are_nonempty_and_flag_theta_qcd_and_lattice():
    res = pc.ed_honest_residuals()
    assert len(res) >= 3
    joined = " ".join(res).lower()
    assert "theta_qcd" in joined
    assert "lattice" in joined           # the load-bearing caveat is stated
    assert "colour" in joined or "color" in joined
