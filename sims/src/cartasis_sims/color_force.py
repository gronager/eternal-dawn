r"""Effective forces from soliton overlaps: the colour channels and the residual force.

The vision: the particles are soliton wavefunctions; their overlaps are effective
potentials; the forces are the channels of those overlaps. This module makes that
concrete for the strong sector.

TWO PIECES, with clearly different status.

(1) The COLOUR-CHANNEL structure -- RIGOROUS, and FORCED by the framework. The soliton
carries the 3-valued antisymmetric label forced by Pauli (the Delta++ argument), i.e.
the SU(3) fundamental. The two-body colour factor <T1.T2> = (C_R - C_1 - C_2)/2 is then
fixed by the SU(3) Casimirs -- identical to QCD, not chosen:

    q-qbar:  singlet 1 -> -4/3 (attractive: MESONS),   octet 8 -> +1/6 (repulsive)
    q-q   :  antitriplet 3bar -> -2/3 (attractive but NOT a singlet: confined diquark),
             sextet 6 -> +1/3 (repulsive)

So free states must be colour SINGLETS: q-qbar (meson, 2-body) or qqq (baryon, 3-body,
the totally antisymmetric singlet) -- and there is NO free 2-quark state, exactly as
observed. The channel structure (which combinations bind) is a consequence of the
forced label.

(2) The radial SHAPE -- MODELLED. The colour-force potential's shape (Coulomb-like at
short range, confining linear at long range -- the Cornell form) comes from the
non-abelian connection dynamics, whose string tension is genuinely lattice-scale and
NOT computed here. We use a Cornell shape as illustration; only the colour FACTORS
multiplying it are rigorous.

(3) The RESIDUAL force between colour-singlet composites -- the "force from overlaps".
Two singlets feel a leftover force from the overlap of their fields (the van der Waals
of the colour force / sigma-exchange between solitons): an attractive Yukawa set by the
condensate field (mass m_sigma). This is the nuclear-force analogue, computable from
the soliton's own field, and it is what the overlap picture delivers directly.
"""

from __future__ import annotations

import numpy as np

# SU(3) quadratic Casimirs
_C2 = {"1": 0.0, "3": 4.0 / 3, "3bar": 4.0 / 3, "8": 3.0, "6": 10.0 / 3, "6bar": 10.0 / 3}

# two-body channels: (combined rep, rep1, rep2, label)
CHANNELS = {
    "qqbar_singlet": ("1", "3", "3bar", "q-qbar singlet (meson)"),
    "qqbar_octet": ("8", "3", "3bar", "q-qbar octet"),
    "qq_antitriplet": ("3bar", "3", "3", "q-q antitriplet (diquark)"),
    "qq_sextet": ("6", "3", "3", "q-q sextet"),
}


def color_factor(channel):
    """<T1.T2> = (C_R - C_1 - C_2)/2 for a two-body colour channel (negative=attractive).
    Forced by the SU(3) fundamental label -- identical to QCD."""
    R, r1, r2, _ = CHANNELS[channel]
    return (_C2[R] - _C2[r1] - _C2[r2]) / 2.0


def cornell(r, alpha=0.4, sigma=1.0):
    """Modelled colour-force shape: -alpha/r + sigma r (Coulomb + linear confinement).
    Only illustrative -- the string tension sigma is lattice-scale, not computed here."""
    r = np.asarray(r, dtype=float)
    return -alpha / r + sigma * r


def channel_potential(r, channel, alpha=0.4, sigma=1.0):
    """V_channel(r) = <T1.T2> * Cornell(r): the colour factor (rigorous) times the
    modelled shape. Attractive (binding) iff the colour factor is negative."""
    return color_factor(channel) * cornell(r, alpha, sigma)


def binds(channel):
    """Does the channel attract (negative colour factor)?"""
    return color_factor(channel) < 0


def residual_yukawa(r, g=1.0, m_sigma=1.0):
    """The residual force between two colour SINGLETS from the overlap of their
    condensate (sigma) fields -- an attractive Yukawa, the nuclear-force analogue.
    This is the 'force from overlaps' the soliton picture delivers directly."""
    r = np.asarray(r, dtype=float)
    return -(g**2) * np.exp(-m_sigma * r) / (4 * np.pi * r)


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
