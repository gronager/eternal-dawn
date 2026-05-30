r"""Black-hole / cosmology identities for the "universe = black-hole interior" thesis.

These are the cheapest, most decisive "Phase 0" checks: take the universe's
*measured* numbers and ask what they imply if we sit inside a (rotating) black
hole. Nothing here assumes the Cartasis mechanism; it only tests consistency.

Key results (derivations in the module functions' docstrings):

* Schwarzschild radius of the Hubble-sphere mass equals the Hubble radius
  exactly when the universe is spatially flat:  R_s / R_H = Omega.
* The maximal causal spin parameter of that mass is of order unity:
  a*_max = c^2 R_H / (G M) = 2 / Omega ~ 2.
* The Hawking temperature of the Hubble-sphere mass is ~10^-30 K, thirty orders
  of magnitude below the CMB -- so "CMB = Hawking" cannot mean present-horizon
  radiation; it must be parent-horizon radiation from the bounce epoch.
"""

from __future__ import annotations

import math

from . import constants as k


# ---------------------------------------------------------------------------
# Geometry / mass
# ---------------------------------------------------------------------------
def schwarzschild_radius(M: float) -> float:
    r"""R_s = 2 G M / c^2."""
    return 2.0 * k.G * M / k.c**2


def hubble_radius(H: float = k.H0) -> float:
    r"""R_H = c / H."""
    return k.c / H


def critical_density(H: float = k.H0) -> float:
    r"""rho_crit = 3 H^2 / (8 pi G)."""
    return 3.0 * H**2 / (8.0 * math.pi * k.G)


def mass_within_hubble(H: float = k.H0, Omega: float = 1.0) -> float:
    r"""Mass-energy inside the Hubble sphere at density rho = Omega * rho_crit.

    M = (4/3) pi R_H^3 (Omega rho_crit).
    """
    R_H = hubble_radius(H)
    return (4.0 / 3.0) * math.pi * R_H**3 * Omega * critical_density(H)


def rs_over_rh(H: float = k.H0, Omega: float = 1.0) -> float:
    r"""R_s / R_H for the Hubble-sphere mass.

    With M = (4/3) pi R_H^3 Omega rho_crit and rho_crit = 3 H^2 / 8 pi G,

        R_s = 2 G M / c^2 = 2 G / c^2 * (4/3) pi (c/H)^3 * Omega * 3 H^2/(8 pi G)
            = Omega * (c / H) = Omega * R_H.

    So a spatially flat universe (Omega = 1) sits exactly at its own
    Schwarzschild radius -- a black-hole property, read straight off the
    measured flatness. Returns the computed ratio (should equal Omega).
    """
    return schwarzschild_radius(mass_within_hubble(H, Omega)) / hubble_radius(H)


# ---------------------------------------------------------------------------
# Spin
# ---------------------------------------------------------------------------
def spin_parameter(J: float, M: float) -> float:
    r"""Dimensionless Kerr spin a* = c J / (G M^2).  Sub-extremal iff a* < 1."""
    return k.c * J / (k.G * M**2)


def max_spin_parameter(H: float = k.H0, Omega: float = 1.0) -> float:
    r"""Maximal *causal* spin parameter of the Hubble-sphere mass.

    Take the largest angular momentum compatible with causality, J ~ M c R_H
    (all the mass at the horizon moving at ~c). Then

        a*_max = c J / (G M^2) = c^2 R_H / (G M) = 2 / Omega,

    using G M / c^2 = R_s / 2 = (Omega R_H)/2. For Omega = 1, a*_max = 2:
    a modest net rotation leaves the universe comfortably sub-extremal (a* < 1),
    consistent with a genuine, cosmic-censored Kerr interior.
    """
    M = mass_within_hubble(H, Omega)
    R_H = hubble_radius(H)
    J_max = M * k.c * R_H
    return spin_parameter(J_max, M)


# ---------------------------------------------------------------------------
# Temperatures
# ---------------------------------------------------------------------------
def hawking_temperature(M: float) -> float:
    r"""T_H = hbar c^3 / (8 pi G M kB)."""
    return k.hbar * k.c**3 / (8.0 * math.pi * k.G * M * k.kB)


def gibbons_hawking_temperature(H: float = k.H0) -> float:
    r"""de Sitter / cosmological-horizon temperature  T_GH = hbar H / (2 pi kB).

    Note the identity (when R_s = R_H): T_H(M_Hubble) = T_GH / 2.
    """
    return k.hbar * H / (2.0 * math.pi * k.kB)


# ---------------------------------------------------------------------------
# Membrane flip/filter (algebraic, 0-D)
# ---------------------------------------------------------------------------
def filter_fraction(omega_b: float = k.Omega_b_hsq,
                    omega_c: float = k.Omega_c_hsq) -> float:
    r"""Filter pass-fraction f if dark matter is filter-rejected counterpart matter.

    Dark-matter-to-baryon ratio = (1 - f) / f = omega_c / omega_b, so

        f = omega_b / (omega_b + omega_c).

    Planck gives omega_c/omega_b ~ 5.4, hence f ~ 0.16 ~ 1/6: the membrane would
    have to pass ~1/6 of the crossing matter and reject ~5/6.
    """
    return omega_b / (omega_b + omega_c)
