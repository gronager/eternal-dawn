r"""Schwarzschild--de Sitter thermodynamics: why an OGU evaporates, and what happens
to its interior when the membrane goes.

This corrects a sloppy earlier claim ("the OGU reaches thermal equilibrium with the
void, so growth stops") and resolves the puzzle it created ("if it's in equilibrium,
why does it evaporate?"). The honest treatment uses the TWO horizons of
Schwarzschild--de Sitter (SdS).

SdS geometry (units G = c = H_Lambda = 1): f(r) = 1 - 2 mu/r - r^2, with mu the mass.
For 0 < mu < mu_Nariai = 1/(3 sqrt 3) ~ 0.192 there are two horizons: a black-hole
horizon r_b and a cosmological horizon r_c > r_b, with surface-gravity temperatures
T = |f'(r)|/4pi. The crucial fact:

    T_b > T_c   ALWAYS   (equal only at the Nariai point, where r_b = r_c = 1/sqrt3).

So the black hole is always HOTTER than its own cosmological horizon. There is never
true equilibrium (only the measure-zero, unstable Nariai coincidence). Net radiation
therefore always flows from the BH to the cosmological horizon: the OGU slowly
EVAPORATES while the cosmological horizon grows. Near Nariai the temperature gap
T_b - T_c -> 0, so the evaporation becomes arbitrarily slow -- a "lukewarm",
near-degenerate plateau that mimics equilibrium but is not it. That is why the OGU
both (a) caps at ~the Nariai mass (the geometric maximum -- no larger BH fits in this
dS) and (b) nonetheless evaporates: the cap is a GEOMETRIC bound, not a thermal
balance, and the residual T_b > T_c drains it over ~1e142 s.

The interior, and "does time inside stop?" No -- and "no discontinuities" is the
guide. The membrane is the baby's PAST cosmological horizon (its Big-Bang surface),
not its future; the baby's own future runs to its own de Sitter infinity, eternally
(recursion.py). Exterior evaporation and interior eternity are duals across the
horizon, not a contradiction: the horizon is a null surface, and an infinite interior
proper time maps to the finite exterior horizon lifetime. By the (exterior) time the
parent's horizon has evaporated, the baby interior has -- on its own clock -- already
expanded to its de Sitter future infinity, i.e. BECOME a flat, heat-dead void hosting
its own sub-OGUs. So the parent's evaporation is the exterior face of the interior
reaching its natural conformal infinity (the flat void); nothing's clock is stopped,
no nested structure is destroyed. The inner OGUs were always inside the baby's own
void = its far future; when the parent horizon closes, that void is simply a
free-standing region of the eternal substrate, carrying its sub-structure with it.
The umbilical channel pinches off (the connection ends); the universe does not.
"""

from __future__ import annotations

import numpy as np

MU_NARIAI = 1.0 / (3.0 * np.sqrt(3.0))      # max mass in units c^3/(G H_Lambda)
R_NARIAI = 1.0 / np.sqrt(3.0)               # degenerate horizon radius


def sds_horizons(mu: float):
    """Black-hole and cosmological horizon radii (r_b, r_c) of Schwarzschild-de
    Sitter for mass mu (units c=G=H=1), roots of r^3 - r + 2 mu = 0. Returns
    (r_b, r_c) with r_b <= r_c, or None above the Nariai mass."""
    roots = np.roots([1.0, 0.0, -1.0, 2.0 * mu])
    real = sorted(r.real for r in roots if abs(r.imag) < 1e-9 and r.real > 0)
    if len(real) < 2:
        return None
    return real[0], real[-1]


def _surface_temperature(mu: float, r: float) -> float:
    """T = |f'(r)|/(4 pi), f'(r) = 2 mu/r^2 - 2 r (units H = 1, so T in units H/2... )."""
    return abs(2.0 * mu / r**2 - 2.0 * r) / (4.0 * np.pi)


def sds_temperatures(mu: float):
    """(T_b, T_c): black-hole-horizon and cosmological-horizon temperatures."""
    h = sds_horizons(mu)
    if h is None:
        return None
    rb, rc = h
    return _surface_temperature(mu, rb), _surface_temperature(mu, rc)


def temperature_ratio(mu: float) -> float:
    """T_b/T_c >= 1: the black hole is always hotter than its cosmological horizon
    (=1 only at Nariai). The driver of net outflow -> evaporation."""
    t = sds_temperatures(mu)
    if t is None:
        return 1.0
    return t[0] / t[1]


def in_true_equilibrium(mu: float, tol: float = 1e-6) -> bool:
    """True only at the (unstable, measure-zero) Nariai point where T_b = T_c."""
    return abs(temperature_ratio(mu) - 1.0) < tol


def net_outflow(mu: float) -> bool:
    """The BH loses mass whenever T_b > T_c, i.e. for every sub-Nariai mass."""
    return temperature_ratio(mu) > 1.0


def nariai_fraction(M_kg: float, M_nariai_kg: float) -> float:
    """Dimensionless mu/mu_Nariai = M/M_Nariai for a physical mass."""
    return M_kg / M_nariai_kg
