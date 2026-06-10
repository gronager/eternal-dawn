r"""The Horizon, step 1: the full fermion spectrum condensing out of ONE chiral transition.

The genesis modules show the condensate FORM; this shows the SPECTRUM it carries -- the framework's
central claim made runnable. The mass law is

    m(tower, gen)  =  Lambda * c_T * |O_n(s_T)|^2 ,

ONE condensate scale Lambda, a per-tower coupling c_T (the charge/colour grip on the condensate:
coloured > charged lepton, the neutrino weakest), and a per-generation radial overlap O_n(s_T) (the
bag sharpness, the three rungs). So as the universe cools and the condensate v(T) switches on, ALL
twelve fermions acquire mass AT ONE TRANSITION -- the unification, unlike the Standard Model's staged
electroweak-then-QCD -- fanning out into the observed hierarchy. This module drives that emergence:
v(T) rises through the transition (the same S-curve the genesis quench produces), and every fermion
mass m(T) = v(T) * m_rel switches on together, spread across c_T (the towers) and s_T (the rungs).

HONEST: this is the mass STRUCTURE made into a cascade -- RELATIVE masses (absolute MeV needs Lambda
= g, b0), and the inter-tower factors + the per-tower s_T are still the lattice residual (we use the
observed spread to draw the fan; the framework's content is that the spread FACTORISES as c_T x rung).
It is the first Horizon step: the Standard Model condensing out of the bounce, parametrised."""
from __future__ import annotations

import numpy as np

from . import fermion_masses as fm

# the four towers, the colour/charge that sets c_T, and the observed masses (relative -- absolute MeV
# is owed to Lambda). Neutrinos illustrative (only Delta m^2 measured).
TOWERS = ("up-quark", "down-quark", "charged-lepton", "neutrino")
_CHARGE = {  # (n_colour, |Q|) -- the gross "grip on the condensate"; sets the c_T ORDERING
    "up-quark": (3, 2 / 3), "down-quark": (3, 1 / 3), "charged-lepton": (1, 1), "neutrino": (1, 0),
}


def relative_masses():
    """The 12 fermion masses as the framework factorises them: m = c_T * rung_n, normalised so the
    electron = 1. Returns {tower: {'m': [3], 'c_T': float, 's_T': float, 'rungs': [3]}}. The c_T (tower
    scale) is the tower's geometric-mean mass; the rungs are the within-tower ladder (the generation
    hierarchy). Drawn from the observed spread -- the framework's claim is the FACTORISATION."""
    e = fm.OBSERVED["charged-lepton"][0]
    out = {}
    for t in TOWERS:
        m = np.array(fm.OBSERVED[t], dtype=float) / e          # relative to the electron
        c_T = float(np.exp(np.mean(np.log(m))))                # tower scale = geometric mean
        rungs = m / c_T                                        # within-tower ladder (generations)
        s_T = fm.fit_source_size(t)[0] if t != "neutrino" else float("nan")
        out[t] = {"m": m, "c_T": c_T, "s_T": s_T, "rungs": rungs}
    return out


def condensate_curve(temps, Tc=0.5, width=0.06):
    """The chiral order parameter v(T) switching on through the transition: a tanh S-curve (the shape
    the genesis quench's mean order produces), 0 in the hot symmetric phase, 1 in the cold vacuum."""
    temps = np.asarray(temps, dtype=float)
    return 0.5 * (1.0 + np.tanh((Tc - temps) / width))


def emergence(n_T=200, T_hot=1.0, T_cold=0.0, Tc=0.5, width=0.06):
    """The spectrum switching on as the universe cools: returns (temps, v, masses) where masses is
    {tower: (n_T, 3)} = v(T) * m_rel -- every fermion acquiring mass AT ONE transition, fanning into
    the hierarchy. Cooling runs T_hot -> T_cold (left to right in the figure)."""
    temps = np.linspace(T_hot, T_cold, n_T)
    v = condensate_curve(temps, Tc=Tc, width=width)
    rel = relative_masses()
    masses = {t: v[:, None] * rel[t]["m"][None, :] for t in TOWERS}
    return temps, v, masses, rel
