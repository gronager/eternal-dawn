"""Horizon 3d: the unified genesis -- all matter from one transition, with the lattice re-run hook."""
import numpy as np
from cartasis_sims import genesis_horizon as gh


def test_full_matter_has_all_sectors():
    m = gh.full_matter()
    assert len(m) == 14                                    # 8 baryons + 3 charged leptons + 3 neutrinos
    sectors = {s for _, s in m.values()}
    assert sectors == {"octet", "decuplet", "lepton", "neutrino"}
    assert m["N"][0] > m["e"][0] > m["nu1"][0]             # the matter hierarchy spans the spectrum


def test_genesis_switches_all_on_at_one_transition():
    g = gh.genesis(n_T=200)
    assert g["v"][0] < 0.02 and g["v"][-1] > 0.98          # condensate off (hot) -> on (cold)
    for sp, (mfin, _) in g["matter"].items():
        assert g["masses"][sp][0] < 0.05 * mfin + 1e-30    # massless hot
        assert np.isclose(g["masses"][sp][-1], mfin)       # full mass cold


def test_lattice_rerun_hook():
    # the owed inputs: scale (Lambda) anchors the nucleon; s_T re-rungs the leptons
    g = gh.genesis(scale_MeV=939.0, s_T=0.45)
    assert np.isclose(g["matter"]["N"][0], 939.0, rtol=1e-6)     # scale anchored on N
    base = gh.genesis()["matter"]["mu"][0]
    assert g["matter"]["mu"][0] != base                          # a different s_T moves the generations
    assert np.isclose(g["matter"]["tau"][0], 1776.86, rtol=1e-6) # tau stays anchored
