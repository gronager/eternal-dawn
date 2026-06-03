"""Tests for the no-fit charge/colour structure that distinguishes the fermions."""

from cartasis_sims import ab_initio_charges as ai


def test_neutrino_is_forced_massless_no_fit():
    # a total-singlet (RH) neutrino has zero condensate handle -- light by its charges
    assert ai.neutrino_is_forced_massless() == 0.0


def test_ordering_is_quarks_then_lepton_then_neutrino():
    order = ai.tower_handles()["ordering"]
    assert order[0] in ("up-quark", "down-quark")
    assert order[1] in ("up-quark", "down-quark")
    assert order[2] == "charged-lepton"
    assert order[3] == "neutrino-RH"


def test_quarks_sit_above_leptons_via_colour():
    H = ai.tower_handles()["handles"]
    assert H["up-quark"] > H["charged-lepton"]
    assert H["down-quark"] > H["charged-lepton"]
    # and the excess is dominated by colour (the handle leptons lack)
    assert ai.colour_contribution_fraction("up-quark") > 0.8
    assert ai.colour_contribution_fraction("charged-lepton") == 0.0


def test_handles_use_only_charges_and_gauge_couplings():
    # changing a charge changes the handle; there is no fitted per-tower knob
    base = ai.condensate_handle("charged-lepton")
    # a hypothetical neutral lepton (Q=0) would couple less -> smaller handle
    import copy
    saved = copy.deepcopy(ai.QUANTUM_NUMBERS["charged-lepton"])
    try:
        ai.QUANTUM_NUMBERS["charged-lepton"]["Q"] = 0.0
        assert ai.condensate_handle("charged-lepton") < base
    finally:
        ai.QUANTUM_NUMBERS["charged-lepton"] = saved


def test_what_is_owed_names_the_lattice_residue():
    w = ai.what_is_owed()
    assert any("neutrino" in s for s in w["ab_initio_no_fit"])
    joined = " ".join(w["owed_to_lattice_L4"]).lower()
    assert "absolute" in joined and "0.511" in joined
