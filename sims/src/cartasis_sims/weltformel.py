r"""The Weltformel spectrum in physical units: all 15 elementary masses, no fitting.

This is the climax assembly. It takes the no-fit pieces of Part~II and combines them with a
SINGLE global scale (not a per-tower anchor):

  * the generation factor O(n) = the substrate overlap of the n-th radial level
    (substrate_overlap.py), one shared well depth, NO fit;
  * the tower factor h(T) = the charge/colour handle (ab_initio_charges.py), pure group theory
    plus the measured gauge couplings, NO fit;
  * ONE global scale Lambda, fixed by the condensate (the heaviest fermion is the unsuppressed
    ground state, m_top ~ v/sqrt2 ~ 173 GeV -- the top Yukawa is ~1), which is the one number
    owed to the lattice (L5 / the transmuted scale of Eq. 11.1).

    m(T, n) = Lambda * h(T)/h_max * O(n)/O_max .

Then W, Z (= (1/2) g v, (1/2) sqrt(g^2+g'^2) v, using the measured gauge couplings) and the
composite Higgs (the walking pseudo-dilaton, ~0.5 v). Fifteen masses, in GeV, from the
equations and the constants -- with exactly ONE overall scale and NO per-particle fitting. The
gross structure is exact (three generations, geometric ladders, quarks above leptons, the
neutrino seesaw-suppressed -- see below); the absolute values are within one to two orders of
magnitude; the bosons land within a percent. The decimals (the inter-tower splittings, the exact
spreads) are the lattice's (Appendix).

NEUTRINOS. A Dirac mass needs BOTH a left- and a right-handed field to grip the condensate. The
right-handed neutrino is a total singlet (colour 0, charge 0, weak 0), so the leading Dirac handle
vanishes and the physical mass is the seesaw value m_nu ~ m_D^2 / M (M a heavy Majorana scale):
tiny, sub-eV, NOT exactly zero, and far below every charged fermion. We carry the leading-order
value as zero and treat the sub-eV seesaw floor (the scale M owed) as the prediction.
"""

from __future__ import annotations

import numpy as np

from . import ab_initio_charges as ai
from . import substrate_overlap as su
from . import electroweak_masses as ew

TOWERS = ["up-quark", "down-quark", "charged-lepton", "neutrino"]
LABELS = {
    "up-quark": ["u", "c", "t"], "down-quark": ["d", "s", "b"],
    "charged-lepton": ["e", "mu", "tau"], "neutrino": ["nu1", "nu2", "nu3"],
}
# PDG masses in GeV (neutrinos: representative ~eV scale, illustrative only)
OBSERVED_GEV = {
    "up-quark": [2.16e-3, 1.273, 173.0],
    "down-quark": [4.7e-3, 93.5e-3, 4.18],
    "charged-lepton": [0.511e-3, 105.66e-3, 1.77693],
    "neutrino": [1e-11, 9e-11, 5e-10],
}
OBSERVED_BOSONS = {"W": 80.37, "Z": 91.19, "H": 125.20}
TOP_GEV = 173.0           # the global scale: heaviest fermion = condensate scale (Yukawa ~1)


def fermion_matrix(top_GeV=TOP_GEV, depth=6.0):
    """The 3x4 fermion mass matrix in GeV, purely from the equations + ONE global scale.
    Generation factor: substrate overlap (no fit). Tower factor: charge handle (no fit).
    The neutrino's leading Dirac mass vanishes (the right-handed partner is a total singlet),
    so its row is zero at leading order; the physical mass is the sub-eV seesaw floor (M owed)."""
    gen = su.substrate_overlap_masses(depth=depth)        # ascending [g1,g2,g3], no fit
    gen = gen / gen.max()
    handles = {t: ai.condensate_handle(t) for t in ("up-quark", "down-quark", "charged-lepton")}
    handles["neutrino"] = ai.condensate_handle("neutrino-RH")   # total singlet -> 0
    hmax = max(handles.values())
    # one global scale so the single largest entry (up gen-3 = top) equals top_GeV
    scale = top_GeV / ((handles["up-quark"] / hmax) * gen[-1])
    return {t: scale * (handles[t] / hmax) * gen for t in TOWERS}


def bosons():
    """W, Z (tree-level, from the condensate scale and the measured gauge couplings) and the
    composite Higgs (walking pseudo-dilaton ~ 0.5 v). In GeV."""
    return {"W": ew.m_W(), "Z": ew.m_Z(), "H": ew.higgs_walking_dilaton()}


# The composite ('techni-hadron') family the picture predicts ABOVE the Higgs. If the 125 GeV
# scalar is a torsiton pair, then -- like QCD's meson tower (sigma, rho, a1, ...) -- there must
# be heavier composites. Rough scale set by the walking bound M_V >~ 13 f_pi (S < 0.1, Ch. 8),
# with v ~ 246 GeV: vectors at a few TeV, axials higher. Estimates only; exact owed to lattice.
COMPOSITE_FAMILY = [
    # (name, symbol, rough mass estimate in GeV, observed?)
    ("scalar (Higgs)", "h_0", ew.higgs_walking_dilaton(), True),
    ("heavier scalar", "h'", 1500.0, False),
    ("vector (techni-rho)", "rho_T", 3000.0, False),
    ("axial-vector (techni-a_1)", "a_T", 5000.0, False),
]


def composite_resonances():
    """The predicted composite (torsiton-pair) family. The Higgs is the lightest member
    (found); the rest are TeV-scale resonances not yet observed -- the sharp collider test
    (the Standard Model predicts a desert above the Higgs; this framework predicts company)."""
    return COMPOSITE_FAMILY


def residual_factors():
    """Per-particle discrepancy max(pred,obs)/min(pred,obs) for the charged fermions and the
    bosons (neutrinos excluded -- predicted zero, no finite ratio)."""
    out = {}
    mat = fermion_matrix()
    for t in ("up-quark", "down-quark", "charged-lepton"):
        for lab, p, o in zip(LABELS[t], mat[t], OBSERVED_GEV[t]):
            out[lab] = max(p, o) / min(p, o)
    bm = bosons()
    for k in ("W", "Z", "H"):
        out[k] = max(bm[k], OBSERVED_BOSONS[k]) / min(bm[k], OBSERVED_BOSONS[k])
    return out


def summary():
    return {"fermions_GeV": fermion_matrix(), "bosons_GeV": bosons(),
            "worst_fermion_factor": max(v for k, v in residual_factors().items()
                                        if k not in ("W", "Z", "H")),
            "worst_boson_factor": max(residual_factors()[k] for k in ("W", "Z", "H"))}
