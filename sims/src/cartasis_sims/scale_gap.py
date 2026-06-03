r"""The scale gap: dimensional transmutation of the spin-gravity coupling.

The Weltformel soliton (ab_initio_spectrum.py) comes out at the Planck mass, ~10^17 too
heavy. This module shows that the gap is NOT a fine-tuning but DIMENSIONAL TRANSMUTATION:
an asymptotically-free coupling, weak at M_Pl, runs to strong coupling and condenses at

    Lambda = M_Pl * exp(-2 pi / (b0 g_T^2)) ,

exponentially below the Planck scale. The torsiton binds at Lambda, not M_Pl.

THE SECOND COUPLING. In pure Einstein-Cartan, torsion is algebraic -- glued to G -- so the
only scale is M_Pl. The natural extension (Poincare gauge gravity) lets the SPIN CONNECTION
PROPAGATE, carrying its own gauge coupling g_T, independent of Newton's G. That coupling is
the 'other coupling that quantizes the spin-gravity field to a reasonable scale': if it is
asymptotically free, its IR condensation scale Lambda is exponentially below M_Pl, and an
ORDINARY value g_T^2 ~ O(0.1) reproduces the observed ~10^17 hierarchy. The same mechanism
already makes Lambda_QCD sit ~10^19 below M_Pl from an O(1/30) coupling -- nothing exotic.

This module computes Lambda(g_T), inverts for the g_T that hits the electroweak scale (and
shows it is unremarkable), runs the one-loop coupling from M_Pl down to Lambda, and rescales
the Weltformel matrix by Lambda/M_Pl to show the absolute numbers fall into the right
ballpark once the scale is transmuted. The remaining structure (generation spread) is the
substrate-overlap refinement, separate from the scale.
"""

from __future__ import annotations

import numpy as np

from . import ab_initio_spectrum as sp

M_PL_GEV = sp.planck_mass_GeV()
EW_SCALE_GEV = 173.0        # the torsiton scale target: top ~ v/sqrt2 (the heaviest fermion)


def transmuted_scale(g2, b0=7.0, M_high=M_PL_GEV):
    """The IR condensation scale from dimensional transmutation:
    Lambda = M_high * exp(-2 pi / (b0 * g2)), g2 = g_T^2 the spin-gravity gauge coupling."""
    return M_high * np.exp(-2.0 * np.pi / (b0 * g2))


def g2_for_scale(Lambda=EW_SCALE_GEV, b0=7.0, M_high=M_PL_GEV):
    """Invert: the coupling g_T^2 needed to transmute M_high down to Lambda. The headline --
    how ORDINARY a number reproduces the 10^17 hierarchy."""
    return 2.0 * np.pi / (b0 * np.log(M_high / Lambda))


def hierarchy_is_ordinary(b0_values=(1.0, 3.0, 7.0, 11.0)):
    """For a range of beta-function coefficients b0, the g_T^2 (and alpha_T = g_T^2/4pi)
    that hits the electroweak scale. All come out O(0.01-0.2): the 10^17 gap is the
    exponential of an unremarkable coupling, not a tuning."""
    out = {}
    for b0 in b0_values:
        g2 = g2_for_scale(b0=b0)
        out[b0] = {"g_T^2": g2, "alpha_T": g2 / (4 * np.pi),
                   "ln(M_Pl/Lambda)": float(np.log(M_PL_GEV / EW_SCALE_GEV))}
    return out


def running_coupling(mu_GeV, g2_at_Lambda=4.0 * np.pi, b0=7.0, Lambda=EW_SCALE_GEV):
    """One-loop running alpha_T(mu) = g_T^2(mu)/4pi of an asymptotically-free spin-gravity
    coupling: strong (~1) at Lambda, weak (~few %) at M_Pl. Anchored at Lambda."""
    # 1/alpha(mu) = 1/alpha(Lambda) + (b0/2pi) ln(mu/Lambda)
    inv_a_L = 4.0 * np.pi / g2_at_Lambda
    inv_a = inv_a_L + (b0 / (2 * np.pi)) * np.log(np.asarray(mu_GeV) / Lambda)
    return 1.0 / inv_a


def rescaled_matrix_GeV(Lambda=EW_SCALE_GEV):
    """The Weltformel 3x4 matrix, rescaled from the Planck scale to the transmuted scale
    Lambda. The heaviest entry (top) is pinned to Lambda; everything else follows the
    ab-initio structure -- now in the observed ballpark, no longer Planck."""
    fm = sp.fermion_matrix_GeV()
    M_Pl = fm["M_Pl_GeV"]
    factor = Lambda / M_Pl
    return {t: fm["matrix_GeV"][t] * factor for t in sp.TOWERS}, fm["M_Pl_GeV"], Lambda


def scale_summary():
    """One-call: the required coupling (ordinary), and the rescaled top/electron to show the
    transmuted matrix lands in the right ballpark."""
    g2 = g2_for_scale()
    mat, M_Pl, Lam = rescaled_matrix_GeV()
    return {
        "M_Pl_GeV": M_Pl,
        "target_scale_GeV": Lam,
        "ln_hierarchy": float(np.log(M_Pl / Lam)),
        "g_T^2_needed": g2,
        "alpha_T_needed": g2 / (4 * np.pi),
        "rescaled_top_GeV": float(mat["up-quark"][-1]),
        "rescaled_electron_GeV": float(mat["charged-lepton"][0]),
        "comment": "an O(0.01-0.1) coupling, dimensionally transmuted, turns M_Pl into the "
                   "electroweak scale -- the 10^17 gap is exp(-40), not a tuning",
    }
