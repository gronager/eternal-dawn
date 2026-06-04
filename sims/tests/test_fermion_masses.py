"""Tests for the best-effort 12-fermion-mass overlap-ladder calculation (L4)."""

import numpy as np

from cartasis_sims import fermion_masses as fm


def test_all_charged_fermions_within_a_small_factor():
    # the headline honest claim: 9 charged masses, 5 orders of magnitude, within ~a few
    assert fm.worst_residual_factor() < 5.0


def test_leptons_are_the_clean_win():
    # the charged-lepton tower comes out best (all within ~30% of observed)
    r = fm.predict_tower("charged-lepton")
    assert np.max(r["residual_factor"]) < 1.5


def test_electron_is_predicted_not_anchored():
    # with the default 'heaviest' anchor, gen-1 (electron) is a PREDICTION; it lands within
    # ~20% of the observed 0.511 MeV -- the tiny electron mass computed, not inserted
    r = fm.predict_tower("charged-lepton", anchor="heaviest")
    assert r["anchor"] == "heaviest"
    assert not np.isclose(r["predicted"][0], r["observed"][0], rtol=1e-6)  # not pinned
    assert r["residual_factor"][0] < 1.25                                  # but close


def test_top_mass_from_condensate_scale():
    # m_top ~ v/sqrt(2): the heaviest fermion is the condensate scale itself (Yukawa ~1),
    # so the absolute scale is derived, not a free anchor
    t = fm.top_from_condensate()
    assert t["off"] < 1.02


def test_heaviest_anchor_is_exact_by_construction():
    # with 'heaviest', gen-3 is the per-tower scale anchor and matches exactly
    for t in ("charged-lepton", "up-quark", "down-quark"):
        r = fm.predict_tower(t, anchor="heaviest")
        assert np.isclose(r["predicted"][-1], r["observed"][-1], rtol=1e-9)


def test_total_span_is_roughly_reproduced():
    # the geometric SPAN of each tower (heaviest/lightest) is reproduced within ~3x
    for t in ("charged-lepton", "up-quark", "down-quark"):
        r = fm.predict_tower(t)
        ratio = r["total_span_pred"] / r["total_span_obs"]
        assert 1 / 3.5 < ratio < 3.5


def test_neutrino_needs_tiny_coupling_not_new_scale():
    nu = fm.neutrino_suppression()
    # lightness is a small coupling ratio (~1e-7), the smallest condensate grip
    assert nu["needed_c_ratio"] < 1e-6
    assert np.isclose(nu["c_nu_over_c_charged"], nu["needed_c_ratio"])


def test_koide_is_reported_and_near_two_thirds():
    Q = fm.koide_Q(fm.OBSERVED["charged-lepton"])
    assert np.isclose(Q, 2 / 3, atol=2e-3)


def test_parameter_count_beats_sm():
    p = fm.n_free_parameters()
    assert p["ed_model"] < p["sm_inserted"]
    assert p["ed_ideal"] == 1
