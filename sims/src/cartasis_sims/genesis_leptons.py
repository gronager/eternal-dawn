r"""Horizon, step 3c: the LEPTONS -- the colourless half of the soliton spectrum.

3b gave the baryons: the topological Skyrmion (B=1), its spectrum organised by FLAVOUR collective
rotation (the octet/decuplet). The leptons are the SAME chiral soliton in a different gauge sector --
COLOURLESS, so unconfined and free (3a) -- and a colourless single-fermion soliton has no flavour
rotor to climb; instead its spectrum is the RADIAL tower of its bound Dirac level. So:

    baryons  =  topological knot  +  flavour rotation   ->  the octet (Lambda, Sigma, Xi ...)
    leptons  =  colourless knot   +  radial excitation  ->  the generations (e, mu, tau)

Same condensate, same configurational mass (m = Lambda c_T |O_n(s_T)|^2), two organising principles
set by colour. The three charged-lepton generations are the three radial rungs (the s_T ladder, from
the soliton well); the neutrinos are the same rungs with the WEAKEST condensate grip -- a colourless,
electrically NEUTRAL knot has the smallest c_T, so the neutrino tower sits far below the rest with no
new scale (the framework's reading of the sub-eV neutrino).

HONEST: the charged-lepton ratios come from the soliton radial ladder + one fitted core size s_T (the
generation hierarchy, still the lattice residual); the neutrino masses are illustrative (only the
Delta m^2 are measured) -- the mechanism (weakest grip -> lightest) is the content. With Lambda from
the lattice the absolute scale follows."""
from __future__ import annotations

import numpy as np

from . import fermion_masses as fm


def charged_lepton_tower():
    """The e, mu, tau generations as the three radial rungs of the colourless soliton (configmass).
    Returns names, predicted/observed masses (MeV), the relative ladder (e=1), and the core size s_T."""
    p = fm.predict_tower("charged-lepton")
    pred = np.asarray(p["predicted"], dtype=float)
    return {"names": ["e", "mu", "tau"], "mass": pred, "obs": np.asarray(p["observed"], dtype=float),
            "rel": pred / pred[0], "s_T": float(p["source_size"])}


def neutrino_tower(c_nu_over_charged=None, m_tau_charged=None):
    """The neutrinos: the SAME radial rungs as the charged leptons, scaled by the WEAKEST condensate
    grip c_nu (a colourless, neutral knot). Mass is linear in the coupling, so the whole tower drops by
    c_nu/c_charged -- the sub-eV scale with no new physics. Returns the (illustrative) masses and the
    suppression. Pass c_nu_over_charged to set the grip; default is the ratio that lands ~0.05 eV."""
    sup = fm.neutrino_suppression()
    ratio = c_nu_over_charged if c_nu_over_charged is not None else sup["c_nu_over_c_charged"]
    cl = charged_lepton_tower()
    masses = cl["rel"] * ratio * (cl["mass"][0] if m_tau_charged is None else m_tau_charged)
    return {"names": ["nu1", "nu2", "nu3"], "mass": np.asarray(masses, dtype=float),
            "c_ratio": float(ratio), "rel": cl["rel"], "reading": sup["reading"]}


def lepton_spectrum():
    """The full lepton sector: the charged-lepton generations (radial rungs) and the neutrino tower
    (weakest grip). The colourless companion to the baryon octet/decuplet of genesis_su3."""
    charged = charged_lepton_tower()
    nu = neutrino_tower()
    return {"charged": charged, "neutrino": nu,
            "organising_principle": "radial rungs (generations), not flavour rotation (octet)",
            "vs_baryons": "same soliton, colourless sector: free (3a), no rotor -> the s_T ladder"}
