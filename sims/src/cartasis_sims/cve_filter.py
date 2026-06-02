r"""Chiral-vortical-effect (CVE) asymmetry estimate at the rotating bounce.

Order-of-magnitude estimate of whether anomalous chiral transport in the rotating
bounce plasma can source a matter--antimatter asymmetry of the observed size.
This is a *scaling* estimate with explicit, clearly-labelled assumptions, not a
derivation; the point is to see whether the mechanism is in the right ballpark or
off by tens of orders.

Physics. At the Cartasis density the matter is a hot relativistic chiral plasma.
A rotating plasma carries, by the chiral vortical effect, an axial current along
its vorticity (rotation axis),

    j5 ~ (T^2/6) * omega        (leading thermal term, natural units),

which is parity-odd and tied to the parent black hole's spin. The controlling
dimensionless small parameter is the vorticity-to-temperature ratio omega/T; the
frozen-in chirality/baryon asymmetry is

    eta ~ C * (omega/T),

with C an O(0.01-10) coefficient absorbing the anomaly coefficient and the
sphaleron conversion efficiency (genuinely uncertain). We compute T and omega at
the bounce and compare eta to the observed baryon-to-photon ratio
eta_obs ~ 6.1e-10.
"""

from __future__ import annotations

import math

from . import constants as k

HBAR_GEV_S = 6.582119569e-25     # hbar in GeV*s
GEV4_TO_JM3 = 2.085e37           # 1 GeV^4 = 2.085e37 J/m^3
GEV_TO_K = 1.160451812e13        # 1 GeV / k_B in kelvin
G_STAR_SM = 106.75               # SM relativistic dof at high T
ETA_OBS = 6.1e-10                # observed baryon-to-photon ratio


def bounce_temperature_GeV(rho_C: float = 1.0e50,
                           g_star: float = G_STAR_SM) -> float:
    """Plasma temperature at rho_C from eps = rho_C c^2 = (g* pi^2/30) T^4."""
    eps_GeV4 = (rho_C * k.c**2) / GEV4_TO_JM3
    return (30.0 * eps_GeV4 / (g_star * math.pi**2)) ** 0.25


def bounce_rate_s(rho_C: float = 1.0e50) -> float:
    """Natural bounce angular rate sqrt(8 pi G rho_C / 3) [s^-1] (= 2 max|H|)."""
    return math.sqrt(8.0 * math.pi * k.G * rho_C / 3.0)


def bounce_rate_GeV(rho_C: float = 1.0e50) -> float:
    return bounce_rate_s(rho_C) * HBAR_GEV_S


def vorticity_over_T(rho_C: float = 1.0e50, g_star: float = G_STAR_SM,
                     spin_fraction: float = 1.0) -> float:
    """omega/T at the bounce; omega = spin_fraction * bounce_rate (<=1 maximal)."""
    return spin_fraction * bounce_rate_GeV(rho_C) / bounce_temperature_GeV(rho_C, g_star)


def asymmetry_estimate(rho_C: float = 1.0e50, g_star: float = G_STAR_SM,
                       spin_fraction: float = 1.0, C: float = 1.0) -> float:
    """eta ~ C * (omega/T)."""
    return C * vorticity_over_T(rho_C, g_star, spin_fraction)


def spinC_needed(rho_C: float = 1.0e50, g_star: float = G_STAR_SM) -> float:
    """The product C * spin_fraction required to match eta_obs."""
    return ETA_OBS / vorticity_over_T(rho_C, g_star, spin_fraction=1.0)


# ----------------------------------------------------------------------------
# Inheritance / dilution track.
#
# The estimate above asks the rotating bounce to MANUFACTURE eta from a
# symmetric bath -- the relevant question only for an original-generation
# universe (OGU), which has no parent to inherit from. A black-hole-born
# universe (BHU_n, n>=1; our case, by the dark-matter+dark-energy anchor) has
# a parent that already completed baryogenesis and sits at its own
# baryon-to-entropy ratio ~ eta_parent. Sphalerons erase B+L but conserve
# B-L, so that net B-L passes through the extrusion as a protected remnant.
#
# Across the bounce baryon number (B-L) is conserved while entropy may be
# multiplied by a dilution factor D >= 1 (irreversible compression/shear, plus
# -- if the CMB is parent Hawking inflow -- new photons in the denominator).
# Since eta = n_B / n_gamma scales as (conserved B-L)/(entropy),
#
#     eta_baby_inherited = eta_parent / D.
#
# The observed asymmetry is then roughly
#
#     eta_baby ~ max( eta_parent / D ,  C * (omega/T) )   [inheritance vs fresh skew]
#
# and D is the crux unknown. It is bounded between two physically meaningful
# anchors computed below: D = 1 (adiabatic/entropy-conserving extrusion ->
# inheritance dominates, sign protected, PURE lineages) and the horizon-
# saturating value D = S_BH / S_radiation (the post-bounce interior thermalizes
# to the parent horizon entropy -> catastrophic dilution, fresh skew forced,
# 41x problem returns for everyone, sign re-rolls -> MIXED lineages).
# ----------------------------------------------------------------------------

def inherited_eta(eta_parent: float = ETA_OBS, dilution: float = 1.0) -> float:
    """eta_baby from inheritance: eta_parent diluted by entropy factor D >= 1."""
    return eta_parent / dilution


def bh_entropy_over_kB(M: float) -> float:
    """Bekenstein-Hawking horizon entropy S_BH/k_B = 4 pi G M^2 / (hbar c)."""
    return 4.0 * math.pi * k.G * M**2 / (k.hbar * k.c)


def radiation_entropy_over_kB(M: float, T_GeV: float) -> float:
    """Entropy S/k_B = (4/3) E/(k_B T) of energy E = M c^2 as thermal radiation
    at the bounce temperature (T given in GeV)."""
    E = M * k.c**2
    kT = T_GeV * GEV_TO_K * k.kB
    return (4.0 / 3.0) * E / kT


def dilution_horizon(M: float, rho_C: float = 1.0e50,
                     g_star: float = G_STAR_SM) -> float:
    """Horizon-saturating dilution D_max = S_BH / S_radiation for a parent of
    mass M, with the radiation entropy evaluated at the bounce temperature.
    Closed form: D_max = 3 pi G M k_B T / (hbar c^3)."""
    T_GeV = bounce_temperature_GeV(rho_C, g_star)
    return bh_entropy_over_kB(M) / radiation_entropy_over_kB(M, T_GeV)


def combined_eta(rho_C: float = 1.0e50, g_star: float = G_STAR_SM,
                 spin_fraction: float = 1.0, C: float = 1.0,
                 eta_parent: float = ETA_OBS, dilution: float = 1.0) -> float:
    """eta_baby ~ max(inherited-diluted, fresh chiral-vortical skew)."""
    return max(inherited_eta(eta_parent, dilution),
               asymmetry_estimate(rho_C, g_star, spin_fraction, C))
