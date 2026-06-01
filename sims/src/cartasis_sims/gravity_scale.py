r"""Gravity-scaled coordinates: the conformal map that makes the foam true to gravity.

The user's construction, made precise. To render objects of wildly different physical
size but comparable GRAVITATIONAL importance at comparable visual size, rescale every
object as if it had a common reference density rho_ref (water, say). A mass M then gets
a gravity-scaled radius

    R_gs(M) = (3 M / (4 pi rho_ref))^(1/3),

and its magnification relative to its true (mean-density rho) size is

    mu = R_gs / R_true = (rho / rho_ref)^(1/3).

The reference density CANCELS in any ratio of two gravity-scaled radii, so the map is
well-defined independent of the choice of rho_ref -- only relative densities matter.

Two levels, two very different outcomes (the key result):

  OGU-in-VOID. An OGU is a horizon-mass black hole, whose MEAN density over its
  Schwarzschild volume is the cosmic critical density ~ the void's dark-energy density
  (the R_s = R_H identity, Phase 0). So rho_OGU / rho_void ~ 1, and the gravity-scaling
  magnification is ~1: an OGU is NOT magnified relative to its void. Gravity-scaled
  coordinates ~ TRUE coordinates at this level, and the dilute verdict survives intact
  -- OGUs remain tiny, round, and ~e^{I/4} horizons apart even after scaling. The
  wallpaper at the OGU level is honestly near-empty whichever coordinates you use.

  BH-inside-a-UNIVERSE. An astrophysical black hole is vastly denser than the cosmic
  mean (rho_BH ~ 1e17 kg/m^3 for a stellar hole), so its gravity-scaling magnification
  is enormous (~1e14). Here the conformal map does real work: it blows up the otherwise
  invisible compact objects so the nested structure (the next generation of universes)
  becomes visible. This is where gravity-scaling earns its name.

So the honest gravity-scaled wallpaper is: at the OGU/void level, essentially the true
dilute scatter (near-empty); at the within-universe level, a richly magnified nest of
compact objects. The two-tier behaviour is itself a prediction of the R_s = R_H
identity that defines a universe.
"""

from __future__ import annotations

import numpy as np

from . import constants as k

RHO_WATER = 1000.0               # kg/m^3, a convenient reference density
RHO_LAMBDA = 6.0e-27             # void dark-energy density
RHO_CRIT = 8.5e-27               # cosmic critical density (z=0)


def gravity_scaled_radius(M, rho_ref: float = RHO_WATER):
    """R_gs = (3 M / (4 pi rho_ref))^(1/3) [m]: the radius M would occupy at rho_ref."""
    return (3.0 * np.asarray(M, float) / (4.0 * np.pi * rho_ref)) ** (1.0 / 3.0)


def schwarzschild_radius(M):
    return 2.0 * k.G * np.asarray(M, float) / k.c**2


def black_hole_mean_density(M):
    """Mean density of a black hole over its Schwarzschild volume, M/((4/3)pi R_s^3).
    Falls as 1/M^2 -- supermassive holes are dilute, stellar holes are dense."""
    Rs = schwarzschild_radius(M)
    return np.asarray(M, float) / ((4.0 / 3.0) * np.pi * Rs**3)


def magnification(rho_true: float, rho_ref: float = RHO_WATER) -> float:
    """Gravity-scaling magnification mu = (rho_true/rho_ref)^(1/3); independent of
    rho_ref in any ratio."""
    return (rho_true / rho_ref) ** (1.0 / 3.0)


def ogu_void_magnification(M_ogu: float = 9.2e52) -> float:
    """The magnification of an OGU relative to its void: (rho_OGU/rho_void)^(1/3) ~ 1,
    because a horizon-mass hole has ~ the cosmic mean density (R_s = R_H)."""
    return (black_hole_mean_density(M_ogu) / RHO_LAMBDA) ** (1.0 / 3.0)


def bh_universe_magnification(M_bh: float, rho_universe: float = RHO_CRIT) -> float:
    """The gravity-scaling magnification of an astrophysical black hole relative to its
    host universe's mean density -- enormous, so nested structure becomes visible."""
    return (black_hole_mean_density(M_bh) / rho_universe) ** (1.0 / 3.0)


def gravity_scaled_separation_ratio(M_a: float, M_b: float) -> float:
    """Ratio of gravity-scaled radii R_gs(M_a)/R_gs(M_b) = (M_a/M_b)^(1/3); the
    rho_ref cancels, so this is what the conformal map actually fixes."""
    return (M_a / M_b) ** (1.0 / 3.0)
