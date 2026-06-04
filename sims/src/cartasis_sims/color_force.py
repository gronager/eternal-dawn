r"""Effective forces from soliton overlaps: the colour channels, and the confinement gap.

The vision: particles are soliton wavefunctions; their overlaps are effective
potentials; the forces are the channels of those overlaps. This module makes that
concrete -- AND is honest about where the overlap picture stops.

(1) The COLOUR-CHANNEL structure -- RIGOROUS, FORCED by the framework. The soliton
carries the Pauli-forced 3-valued antisymmetric label = the SU(3) fundamental, so the
two-body colour factor <T1.T2> = (C_R - C_1 - C_2)/2 is fixed -- identical to QCD:

    q-qbar:  singlet 1 -> -4/3 (attractive: MESONS),   octet 8 -> +1/6 (repulsive)
    q-q   :  antitriplet 3bar -> -2/3 (attractive but NOT a singlet: confined diquark),
             sextet 6 -> +1/3 (repulsive)

Free states must be colour SINGLETS: q-qbar (meson, 2-body) or qqq (baryon, 3-body,
the totally antisymmetric singlet). No free 2-quark state -- exactly as observed.

(2) The one-gluon-exchange (OGE) potential, V_OGE = <T1.T2> * alpha/r, is short-range:
it -> 0 as r -> infinity. So is the residual sigma-exchange between singlets (a Yukawa,
the nuclear-force analogue). **Every overlap/exchange force the soliton picture
produces SCREENS -- it asymptotes to zero, not to a linear ~r.**

(3) THE CONFINEMENT GAP -- owed, and a genuine risk. Real quarks are CONFINED: the
quark-quark potential rises linearly, V ~ sigma r, never going to zero (you cannot
isolate a quark). That linear term is NON-PERTURBATIVE -- a collective flux tube of the
non-abelian connection, NOT a pairwise exchange/overlap. The overlap method here does
NOT produce it; it gives screened (-> 0) forces. So `linear_confinement(r)` below is
shown for CONTRAST as what the framework MUST eventually produce, but it is imposed,
not derived. Whether the gravity-torsion connection actually forms flux tubes (confines,
~r) or merely screens (-> 0) is the deepest open question of the strong sector -- if it
only screens, quarks would not be confined and the picture fails here. It is
lattice-scale, undecided, and honestly flagged.
"""

from __future__ import annotations

import numpy as np

_C2 = {"1": 0.0, "3": 4.0 / 3, "3bar": 4.0 / 3, "8": 3.0, "6": 10.0 / 3, "6bar": 10.0 / 3}

CHANNELS = {
    "qqbar_singlet": ("1", "3", "3bar", "q-qbar singlet (meson)"),
    "qqbar_octet": ("8", "3", "3bar", "q-qbar octet"),
    "qq_antitriplet": ("3bar", "3", "3", "q-q antitriplet (diquark)"),
    "qq_sextet": ("6", "3", "3", "q-q sextet"),
}


def color_factor(channel):
    """<T1.T2> = (C_R - C_1 - C_2)/2 (negative = attractive). Forced by the SU(3)
    fundamental label; identical to QCD."""
    R, r1, r2, _ = CHANNELS[channel]
    return (_C2[R] - _C2[r1] - _C2[r2]) / 2.0


def oge_potential(r, channel, alpha=0.4):
    """One-gluon-exchange (perturbative) potential: V = <T1.T2> * alpha / r.
    Sign is set by the colour factor: singlet attractive, octet repulsive. It is
    SHORT-RANGE -- V -> 0 as r -> infinity (a Coulomb-like exchange)."""
    r = np.asarray(r, dtype=float)
    return color_factor(channel) * alpha / r


def residual_yukawa(r, g=1.0, m_sigma=1.0):
    """Residual force between two colour SINGLETS from the overlap of their condensate
    (sigma) fields -- an attractive Yukawa, the nuclear-force analogue. SCREENED:
    V -> 0 as r -> infinity. This is what the overlap picture delivers directly."""
    r = np.asarray(r, dtype=float)
    return -(g**2) * np.exp(-m_sigma * r) / (4 * np.pi * r)


def linear_confinement(r, sigma=1.0):
    """V = sigma r: the linear confining potential real quarks obey (V -> infinity, never
    zero). Shown for CONTRAST -- it is NON-PERTURBATIVE (a flux tube), NOT produced by
    any overlap/exchange here, and is imposed, not derived. The open question is whether
    the gravity-torsion connection actually generates it (confines) or only screens."""
    r = np.asarray(r, dtype=float)
    return sigma * r


def binds_perturbatively(channel):
    """Does the channel attract at short range (negative colour factor)?"""
    return color_factor(channel) < 0


def asymptotes_to_zero(channel, alpha=0.4):
    """Every overlap/exchange force here screens: |V(r)| -> 0 as r grows."""
    return abs(float(oge_potential(8.0, channel, alpha))) < abs(float(oge_potential(1.0, channel, alpha)))


def casimir_fundamental():
    """Validate <T.T> for the fundamental = C_F = 4/3 via Gell-Mann matrices."""
    l = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
         np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]]),
         np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
         np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
         np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]]),
         np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
         np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]),
         np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3)]
    T = [m / 2 for m in l]
    return float(np.real(sum(t @ t for t in T)[0, 0]))
