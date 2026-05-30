r"""Tier 1b: Oppenheimer--Snyder collapse with a torsion bounce interior.

The 1939 Oppenheimer--Snyder model glues a collapsing homogeneous dust ball
(interior: a closed/flat FRW patch) to an exterior Schwarzschild vacuum across
the stellar surface. In standard GR the interior collapses to a singularity
(a -> 0) behind the horizon. Replace the interior dynamics with the
Einstein--Cartan bounce (Tier 1): the interior now bounces at the Cartasis
density rho_C instead of reaching a singularity. Birkhoff's theorem keeps the
exterior exactly Schwarzschild, so:

* from OUTSIDE the surface freezes at the horizon (Schwarzschild time t -> infinity
  as the areal radius R -> r_s) and its light redshifts to zero -- a black hole;
* from INSIDE the dust reaches rho_C at a tiny areal radius R_min, bounces, and
  re-expands -- an inverse bubble, hidden behind the horizon.

This module computes both sides. The interior bounce reuses `bounce` with w = 0
(pressureless dust). The exterior uses the radial timelike geodesic of the
surface (a shell released from rest at R_0):

    Etilde = sqrt(1 - r_s/R_0),
    (dR/dtau)^2 = r_s/R - r_s/R_0,
    dt/dtau = Etilde / (1 - r_s/R),

so the surface reaches r_s in finite proper time tau but infinite Schwarzschild
time t, and its surface redshift is 1 + z = (1 - r_s/R)^(-1/2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import bounce as bnc
from . import constants as k


@dataclass
class OSSolution:
    # interior (proper time tau, areal radius R in units of R_min = 1)
    tau: np.ndarray
    R_int: np.ndarray
    tau_os: np.ndarray
    R_os: np.ndarray          # GR (no spin): collapses to 0
    r_s_units: float          # horizon areal radius, same units (R_min = 1)
    # exterior view (areal radius in units of r_s)
    R_ext: np.ndarray
    t_ext: np.ndarray         # Schwarzschild time (diverges at horizon)
    z_ext: np.ndarray         # surface redshift 1 + z
    # realistic compactness of the bounce
    realistic: dict = field(default_factory=dict)


def r_min(M: float, rho_C: float = 1.0e50) -> float:
    """Areal radius at the bounce: (4/3) pi rho_C R_min^3 = M."""
    return (3.0 * M / (4.0 * math.pi * rho_C))**(1.0 / 3.0)


def r_schwarzschild(M: float) -> float:
    return 2.0 * k.G * M / k.c**2


def rmin_over_rs(M: float, rho_C: float = 1.0e50) -> float:
    """How deep inside the horizon the bounce sits: R_min / r_s."""
    return r_min(M, rho_C) / r_schwarzschild(M)


def exterior_geodesic(R0_over_rs: float, n: int = 2000, eps: float = 1e-3):
    """Schwarzschild time t(R) and redshift (1+z)(R) for the infalling surface.

    Integrated in proper time tau (no turning-point singularity): for a shell
    released from rest at R0 (units r_s = 1),

        dR/dtau = -sqrt(1/R - 1/R0),   dt/dtau = Etilde / (1 - 1/R),

    stopping at R = 1 + eps. Returns (R, t, 1+z) sampled along the worldline;
    t -> infinity and (1+z) -> infinity as R -> 1 (the frozen star).
    """
    from scipy.integrate import solve_ivp

    R0 = R0_over_rs
    Et = math.sqrt(1.0 - 1.0 / R0)

    def rhs(tau, y):
        R = y[0]
        dR = -math.sqrt(max(1.0 / R - 1.0 / R0, 0.0))
        dt = Et / (1.0 - 1.0 / R)
        return [dR, dt]

    def hit_horizon(tau, y):
        return y[0] - (1.0 + eps)
    hit_horizon.terminal = True
    hit_horizon.direction = -1.0

    sol = solve_ivp(rhs, (0.0, 1.0e6), [R0 * (1.0 - 1e-9), 0.0],
                    events=hit_horizon, rtol=1e-10, atol=1e-13,
                    dense_output=True, max_step=0.05)
    tau_end = sol.t_events[0][0]
    taus = np.linspace(0.0, tau_end, n)
    Y = sol.sol(taus)
    R = Y[0]
    t = Y[1]
    oneplusz = 1.0 / np.sqrt(np.clip(1.0 - 1.0 / R, 1e-300, None))
    return R, t, oneplusz


def simulate_os(a_init: float = 8.0, rs_units: float = 3.0,
                M: float = 9.246e52, rho_C: float = 1.0e50) -> OSSolution:
    """Full Oppenheimer--Snyder bounce: interior (units R_min=1) + exterior.

    `rs_units` places the horizon at R = rs_units in interior units (R_min = 1);
    it is illustrative -- the realistic R_min/r_s is reported in `.realistic`.
    `M`, `rho_C` set the realistic compactness numbers.
    """
    sol = bnc.simulate_bounce(w=0.0, a_init=a_init, t_max=200.0, n_points=4000)
    gr = bnc.simulate_gr_collapse(w=0.0, a_init=a_init, n_points=4000)
    # align GR collapse to start where the interior collapse starts
    tau_os = gr["t"] + sol.t[0]

    R0_over_rs = a_init / rs_units
    R_ext, t_ext, z_ext = exterior_geodesic(R0_over_rs)

    realistic = {
        "M_kg": M, "rho_C": rho_C,
        "r_min_m": r_min(M, rho_C), "r_s_m": r_schwarzschild(M),
        "rmin_over_rs": rmin_over_rs(M, rho_C),
        "rmin_over_rs_stellar": rmin_over_rs(10 * k.M_sun, rho_C),
        "bounce_volume_m3": (4.0 / 3.0) * math.pi * r_min(M, rho_C)**3,
    }
    return OSSolution(tau=sol.t, R_int=sol.a, tau_os=tau_os, R_os=gr["a"],
                      r_s_units=rs_units, R_ext=R_ext, t_ext=t_ext, z_ext=z_ext,
                      realistic=realistic)
