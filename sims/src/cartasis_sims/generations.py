r"""Generations: an arithmetic ladder of internal levels -> a geometric mass hierarchy.

The three generations (e/mu/tau, and the quark towers) are the program's quarantined
gap: no ordinary ladder gives the observed hierarchy, which spans orders of magnitude.
This module shows the mechanism that *could* bridge it -- "Yukawa as overlap."

THE IDEA. Suppose the generations are internal excitation levels of the same soliton
(same colour, charge, weak isospin; different internal rung). The level ENERGIES are
roughly arithmetic (evenly spaced, E_n ~ E_0 + n*dE). But the observable MASS is not
the level energy -- in the soliton picture it is a wavefunction OVERLAP with the
localized mass-giving region (the condensate / the would-be Higgs), and overlaps are
EXPONENTIALLY sensitive to the level. A higher rung sits further out / penetrates a
barrier less, so its overlap with the localized core falls like

    O_n  ~  exp(-c n),      mass_n ~ |O_n|^2 ~ exp(-2 c n).

So an ARITHMETIC ladder in n becomes a GEOMETRIC ladder in mass: ln(mass) linear in n,
each generation a fixed multiple of the next. Equivalently: the "unnaturally tiny"
electron Yukawa (~3e-6) is just exp(-2c*0_offset) -- the exponential of an order-one
overlap suppression, not a fine-tuned small number. That is the qualitative win the
Standard Model lacks: a large hierarchy from O(1) inputs.

WHAT IS AND ISN'T DELIVERED. Computed (here, from the real soliton wavefunctions):
the overlap falls exponentially, producing a geometric hierarchy whose steepness is
set by the source size. Observed: the real fermion masses ARE approximately evenly
spaced in log (geometric-ish), exactly the structure this produces. NOT delivered:
the precise ratios (leptons 207, 17 -- a decreasing, not constant, ratio), which need
the detailed level/overlap structure and are the owed refinement. The mechanism
explains the SHAPE of the hierarchy, not yet its exact numbers.
"""

from __future__ import annotations

import numpy as np

from . import soliton as so

_trapz = getattr(np, "trapezoid", None) or np.trapz

# Standard-Model charged-fermion masses (MeV), PDG-ish, for the log-spacing comparison.
SM_FERMIONS = {
    "leptons (e, mu, tau)": [0.511, 105.66, 1776.86],
    "up-type (u, c, t)": [2.16, 1270.0, 172690.0],
    "down-type (d, s, b)": [4.67, 93.4, 4180.0],
}


def core_overlap(E, source_size=0.25, kind="linear", depth=6.0, width=1.0):
    """Overlap of the level-E internal wavefunction with a localized core source
    chi(r) = exp(-(r/source_size)^2): O = |<chi|G>| / ||psi||."""
    r, G, F = so.wavefunction(E, kind=kind, depth=depth, width=width)
    chi = np.exp(-(r / source_size) ** 2)
    norm = _trapz((G**2 + F**2) * r**2, r)
    o = _trapz(G * chi * r**2, r)
    return abs(o) / np.sqrt(max(norm, 1e-300))


def overlap_masses(levels, source_size=0.25, kind="linear", depth=6.0, width=1.0):
    """mass_n ~ |O_n|^2 for a set of internal levels (the overlap hierarchy)."""
    ov = np.array([core_overlap(E, source_size, kind, depth, width) for E in levels])
    return ov**2


def geometric_factor(masses):
    """The per-generation mass ratio from a log-linear fit: each rung is `factor` x
    the next. A constant factor == a clean geometric hierarchy."""
    n = np.arange(len(masses))
    slope = np.polyfit(n, np.log(masses), 1)[0]
    return float(np.exp(-slope))


def log_spacings(masses):
    """Successive ln(mass) gaps. Equal gaps == geometric. The SM fermions have gaps
    that are equal to within a factor ~2 -- approximately geometric."""
    lm = np.log(np.asarray(masses, dtype=float))
    return np.diff(lm)


def is_approximately_geometric(masses, tol=2.5):
    """True if the ln-mass gaps are equal to within a factor `tol` (geometric-ish)."""
    g = np.abs(log_spacings(masses))
    return float(np.max(g) / np.min(g)) < tol
