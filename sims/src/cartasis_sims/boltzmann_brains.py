r"""Why we are not Boltzmann brains: the three towers and the arrow-priority argument.

Companion to chapters/07-time-arrow.tex, "Why we are not Boltzmann brains".

The naive worry compares three de Sitter nucleation costs (probabilities ~ e^{-I},
I = 2 pi M c^2 / (hbar H_Lambda) = the horizon-entropy deficit):

    brain (~1 kg)        : P ~ 10^(-10^69)
    OG seed (~1e15 kg)   : P ~ 10^(-10^84)
    bare Big Bang        : P ~ 10^(-10^123)   (Penrose's 1-in-10^(10^123))

The brain looks cheapest, so "typical observers are fluctuations". This module
shows, numerically, (a) that the three are double-exponential towers whose inner
exponents are 69 / 84 / 123, (b) that no single-tower rescue (head-count,
recursion, lifetime) can close the gap, and (c) the actual resolution: the brain
and the bare Big Bang are LOW-ENTROPY ORGANIZED states that nucleation cannot
deliver from a stationary, arrowless eigenstate -- only generic high-entropy
blobs (viable OG seeds) can nucleate, and they make their arrow at the bounce.
The interior-brain channel is then killed by SCT's FINITE dark energy: the de
Sitter recurrence time needed to fluctuate one brain (~10^(10^69) yr) hugely
exceeds the dark-energy lifetime (~10^136 yr, parent evaporation), whereas an
eternal-Lambda LambdaCDM future does not have this cutoff.
"""

from __future__ import annotations

import numpy as np

from . import constants as k
from . import void_eigenstate as ve

# Penrose's gravitational-entropy bound for our universe's special initial state.
PENROSE_INNER = 123.0          # bare Big Bang: P ~ 1 in 10^(10^123)
M_BRAIN = 1.4                  # kg
M_SEED = 9.0e14               # kg, minimal self-gravitating OG seed (Ch5)
DE_LIFETIME_YEARS = 1.0e136    # parent-evaporation cutoff of SCT dark energy (~1e53 kg)


def instanton_action(M: float) -> float:
    """de Sitter nucleation action I = 2 pi M c^2 / (hbar H_Lambda) = horizon
    entropy deficit of inserting mass M. P_nucleate ~ e^{-I}."""
    H = ve.hubble_lambda()
    return 2.0 * np.pi * M * k.c**2 / (k.hbar * H)


def tower_inner_exponent(M: float) -> float:
    """The inner exponent x in 'P ~ 1 in 10^(10^x)': I/ln(10) = 10^x."""
    return float(np.log10(instanton_action(M) / np.log(10.0)))


def three_towers() -> dict:
    """Inner exponents of the brain / seed / bare-Big-Bang towers."""
    return {
        "brain": tower_inner_exponent(M_BRAIN),
        "seed": tower_inner_exponent(M_SEED),
        "big_bang": PENROSE_INNER,
    }


def per_particle_mass_vs_config() -> float:
    """Ratio of the de Sitter mass-energy toll to the configuration toll, PER
    nucleon: (2 pi m_N c^2 / hbar H) / (~1 bit). ~5e42 -- mass always dominates
    configuration, so a brain can never be made expensive by its organization."""
    m_N = 1.673e-27                       # nucleon mass, kg
    return float(instanton_action(m_N))   # per-nucleon mass action; config ~ 1


def single_tower_rescues_fail() -> bool:
    """No single-tower correction (head-count ~1e50, recursion ~1e3000, lifetime
    ~1e107) approaches the seed-vs-brain gap of 10^(10^84). Returns True: every
    plausible rescue exponent is < 1e84 (the gap's inner exponent)."""
    head_count = 50.0          # log10 observers per universe
    recursion = 3000.0         # log10 of the whole descendant tree (generous)
    lifetime = 107.0           # log10 (universe lifespan / brain thought) in seconds
    gap_inner = three_towers()["seed"]    # ~84
    return max(head_count, recursion, lifetime) < 10.0**gap_inner


def brain_recurrence_time_years() -> float:
    """de Sitter recurrence time to fluctuate one interior brain, ~e^{I_brain}
    Hubble times, in years. Astronomically larger than any finite dark-energy
    lifetime, so a FINITE-dark-energy universe makes zero interior brains."""
    H = ve.hubble_lambda()
    t_H_years = (1.0 / H) / k.year
    log10_years = instanton_action(M_BRAIN) / np.log(10.0) + np.log10(t_H_years)
    return log10_years            # returns log10(years); the value itself overflows


def finite_dark_energy_beats_brains() -> bool:
    """SCT's dark energy switches off (parent evaporation, ~1e136 yr) eons before
    the brain recurrence time (~10^(10^69) yr), so no interior brain ever forms.
    Eternal-Lambda LambdaCDM lacks this cutoff."""
    log10_brain_time = brain_recurrence_time_years()
    log10_de_life = np.log10(DE_LIFETIME_YEARS)
    return log10_brain_time > 1e60 > log10_de_life
