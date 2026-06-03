r"""Part II/III: the charge/colour structure that distinguishes the fermions -- NO FITTING.

The strict ab-initio question: from the gravity-torsion well, its eigenstates, and the
fermion's CHARGES alone (colour rep, electric charge Q, weak isospin T) -- with no fitted
shape and no fitted per-tower coupling -- do the differences between up-quark, down-quark,
charged lepton, and NEUTRINO fall out? This module separates, honestly, the part that does
(structure, ordering, and the forced masslessness of the neutrino) from the part that does
not (the absolute Yukawa normalisations -- the lattice residue, L4).

THE NO-FIT INPUTS. Only group theory and the measured gauge couplings (which are not fitted
here -- they are the SM's three couplings, which ED holds trace to G via the connection):
  * colour quadratic Casimir C2: triplet (quarks) = 4/3, singlet (leptons) = 0;
  * electric charge Q;
  * weak isospin T (and its Casimir T(T+1)).
A fermion's coupling to the mass-giving condensate is built from THESE, with the gauge
couplings as weights -- no free parameter.

THE FORCED RESULT (clean, no fit). A right-handed neutrino is a total singlet: colour 0,
Q 0, weak 0. It couples to NOTHING that condenses, so its configurational mass is ZERO at
leading order -- the neutrino is light because it carries no handle for the condensate to
grip, forced by its charges. Quarks (coloured) couple to the colour condensate that leptons
cannot, so quarks sit above leptons. These orderings are ab initio.

THE WALL (owed, L4). The ABSOLUTE scale of each coupling to the strongly-coupled condensate
-- hence the actual MeV values and the inter-tower factors (top/bottom ~41, top/tau ~97) --
is a non-perturbative normalisation: the lattice. No laptop computation yields 0.511 MeV
from first principles; that is the open problem. This module therefore reports STRUCTURE and
RELATIVE orderings with no fit, and marks the absolute numbers as owed.
"""

from __future__ import annotations

import numpy as np

# Measured gauge couplings (inputs, NOT fitted): strong, electromagnetic, weak.
ALPHA_S = 0.118        # alpha_s(M_Z)
ALPHA_EM = 1 / 128.0   # alpha_em(M_Z)
ALPHA_W = 0.0337       # alpha_2 = alpha_em / sin^2 theta_W

# Colour quadratic Casimir C2(R): fundamental (triplet) 4/3, singlet 0.
C2_TRIPLET = 4.0 / 3.0
C2_SINGLET = 0.0

# Quantum numbers per fermion tower (left-handed; the RH neutrino is a total singlet).
QUANTUM_NUMBERS = {
    "up-quark":       dict(C2=C2_TRIPLET, Q=+2.0 / 3.0, T=0.5),
    "down-quark":     dict(C2=C2_TRIPLET, Q=-1.0 / 3.0, T=0.5),
    "charged-lepton": dict(C2=C2_SINGLET, Q=-1.0,       T=0.5),
    "neutrino":       dict(C2=C2_SINGLET, Q=0.0,        T=0.5),   # LH; RH partner is singlet
    "neutrino-RH":    dict(C2=C2_SINGLET, Q=0.0,        T=0.0),   # total singlet
}


def condensate_handle(tower, alpha_s=ALPHA_S, alpha_em=ALPHA_EM, alpha_w=ALPHA_W):
    """The fermion's total coupling 'handle' on the condensate, from CHARGES ONLY (no fit):

        H = alpha_s * C2_colour  +  alpha_em * Q^2  +  alpha_w * T(T+1).

    Each term is a gauge interaction the condensate can act through; the weights are the
    measured gauge couplings (inputs, not fits). A total singlet has H = 0 -> no mass."""
    qn = QUANTUM_NUMBERS[tower]
    return (alpha_s * qn["C2"]
            + alpha_em * qn["Q"] ** 2
            + alpha_w * qn["T"] * (qn["T"] + 1.0))


def neutrino_is_forced_massless():
    """The clean ab-initio result: a right-handed (total-singlet) neutrino has H = 0 exactly
    -- no charge couples to the condensate, so its leading configurational mass vanishes.
    Light by its charges, not by a tuned small number. Returns the handle (= 0)."""
    return condensate_handle("neutrino-RH")


def tower_handles():
    """The no-fit coupling handles for the four towers, and their RELATIVE ordering -- the
    structure the charges force, with the absolute normalisation owed to lattice (L4)."""
    towers = ["up-quark", "down-quark", "charged-lepton", "neutrino-RH"]
    H = {t: condensate_handle(t) for t in towers}
    order = sorted(H, key=H.get, reverse=True)
    Hmax = max(H.values())
    return {
        "handles": H,
        "relative": {t: H[t] / Hmax for t in towers},   # normalised, no fit
        "ordering": order,                               # quarks > lepton > neutrino(=0)
        "neutrino_handle": H["neutrino-RH"],
    }


def colour_contribution_fraction(tower):
    """Fraction of a quark's handle that comes from COLOUR (the QCD condensate leptons lack)
    -- the no-fit reason quarks sit above leptons."""
    qn = QUANTUM_NUMBERS[tower]
    colour = ALPHA_S * qn["C2"]
    return colour / condensate_handle(tower) if condensate_handle(tower) > 0 else 0.0


def what_is_owed():
    """The explicit wall: what the charges + eigenstates give ab initio, and what is the
    lattice residue (L4)."""
    return {
        "ab_initio_no_fit": [
            "the generation ladder (eigenstates of the self-consistent well; no shape fit)",
            "the neutrino is light/massless (total singlet -> zero condensate handle; forced)",
            "quarks heavier than leptons (colour handle leptons lack)",
        ],
        "owed_to_lattice_L4": [
            "the ABSOLUTE coupling normalisations -> the actual MeV values (e.g. 0.511 MeV)",
            "the inter-tower factors top/bottom ~41, top/tau ~97 (strong-coupling dynamics)",
            "the up/down isospin splitting and Koide's 2/3",
        ],
        "verdict": "charges + eigenstates fix the STRUCTURE and ordering with no fit (the "
                   "neutrino cleanly); the absolute Yukawa normalisations are the genuine "
                   "non-perturbative residue -- no laptop yields them from first principles",
    }
