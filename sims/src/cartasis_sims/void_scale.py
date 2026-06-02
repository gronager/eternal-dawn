r"""The weight of the void and the de Sitter cap on OGU growth (dilute regime).

Two consequences of the dilute verdict (birth_rate.py: I >> I_crit), both set by the
dark energy of the void.

1. Settlers, not city dwellers (the growth cap). In the dilute regime OGUs never
   impinge on neighbours, so the packed-foam "running out of shared void" does not
   apply. Instead a lone OGU's growth is choked by the de Sitter EXPANSION of its own
   void: it can only gather what lies within the cosmological horizon R_dS = c/H_Lambda,
   and the accelerating expansion carries everything else beyond reach. The rigorous
   statement is the Nariai bound -- the largest black hole that fits in de Sitter has

       M_Nariai = c^3 / (3 sqrt(3) G H_Lambda) ~ 4e52 kg,

   so the OGU mass is CAPPED at ~the Nariai mass ~ a Hubble mass ~ our universe, set by
   dark energy. So a settler still "runs out of nothing": not because a neighbour ate
   it, but because the frontier (the horizon) recedes faster than it can settle it.

2. The void weighs something (the gravity-scaled wallpaper). The vacuum has a nonzero
   zero-point / dark energy density rho_Lambda = Omega_Lambda * 3 H_0^2/(8 pi G) ~
   6e-27 kg/m^3. Nonzero weight means nonzero gravity-scaled volume, so universes do
   NOT touch -- they are round and isolated with real gaps (not polygons). And because
   the regime is dilute (OGUs separated by ~e^{I/4} horizons), the void's gravitational
   content overwhelmingly dominates: true to scale, the supraverse is almost entirely
   void, with vanishingly sparse round island universes. The void weight is what gives
   them room to be round; the dilution makes that room enormous.
"""

from __future__ import annotations

import numpy as np

from . import constants as k

GEV_PER_J = 6.242e9            # GeV per joule


def hubble_lambda() -> float:
    """de Sitter rate from dark energy, H_0 sqrt(Omega_Lambda) [s^-1]."""
    return k.H0 * np.sqrt(k.Omega_Lambda)


def critical_density() -> float:
    """rho_crit = 3 H_0^2/(8 pi G) [kg/m^3]."""
    return 3.0 * k.H0**2 / (8.0 * np.pi * k.G)


def vacuum_density() -> float:
    """The void's weight: rho_Lambda = Omega_Lambda rho_crit [kg/m^3]. Nonzero, so the
    void occupies nonzero gravity-scaled volume -- universes are round, not polygons."""
    return k.Omega_Lambda * critical_density()


def vacuum_energy_gev_m3() -> float:
    """Void energy density in GeV/m^3 (a few proton masses per cubic metre)."""
    return vacuum_density() * k.c**2 * GEV_PER_J


def de_sitter_horizon_mass() -> float:
    """Mass within a de Sitter horizon, M_dS = c^3/(2 G H_Lambda) [kg]."""
    return k.c**3 / (2.0 * k.G * hubble_lambda())


def nariai_mass() -> float:
    """Largest black hole that fits in de Sitter (the OGU growth cap),
    M_Nariai = c^3/(3 sqrt(3) G H_Lambda) [kg] ~ a Hubble mass."""
    return k.c**3 / (3.0 * np.sqrt(3.0) * k.G * hubble_lambda())


def gravity_scaled_void_fraction(ogu_separation_horizons: float) -> float:
    """Fraction of gravitational content in the void rather than in universes, for an
    OGU spacing (in de Sitter horizon radii). Each OGU ~ a horizon mass; the void in
    its domain ~ rho_Lambda * (separation)^3 ~ M_dS * sep^3, so void dominates for
    sep >> 1 (the dilute regime). Returns void/(void+universes) in [0,1)."""
    sep3 = ogu_separation_horizons ** 3
    return sep3 / (sep3 + 1.0)
