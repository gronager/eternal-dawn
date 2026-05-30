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
