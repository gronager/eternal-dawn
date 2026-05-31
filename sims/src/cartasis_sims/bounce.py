r"""Tier 1 (minimal): the Einstein--Cartan torsion bounce.

A homogeneous, isotropic spin fluid in Einstein--Cartan obeys an effective
Friedmann equation in which the averaged fermion spin contributes a negative
quadratic term, so the bare density's gravitational pull is cancelled at a
critical (Cartasis) density:

    H^2 = (8 pi G / 3) rho (1 - rho / rho_C),         (dimensionless: rho_C = 1)

with the continuity equation rho ~ a^{-n}, n = 3(1 + w). At rho = rho_C the
expansion rate vanishes and reverses: a bounce, not a singularity. This is the
leading-order form shared by Poplawski's torsion bounce and loop quantum
cosmology; here we simply integrate it and exhibit the bounce.

We work in dimensionless units (8 pi G / 3 = 1, rho_C = 1, a_min = 1) and
integrate the *second-order* (Raychaudhuri) form, which passes smoothly through
H = 0:

    a'' = a^{1-n}(1 - n/2) + (n - 1) a^{1-2n}.

`physical_timescale` converts the dimensionless time back to seconds for a given
Cartasis density.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from . import constants as k


@dataclass
class BounceSolution:
    t: np.ndarray          # dimensionless time (bounce near t = 0 region)
    a: np.ndarray          # scale factor (a_min = 1 at the bounce)
    H: np.ndarray          # Hubble rate  a'/a
    rho: np.ndarray        # density in units of rho_C
    a_min: float
    rho_max: float         # in units of rho_C (should be ~1)
    w: float


def _n(w: float) -> float:
    return 3.0 * (1.0 + w)


def _accel(a: float, n: float) -> float:
    """a'' from the effective Friedmann + continuity equations."""
    return a**(1.0 - n) * (1.0 - n / 2.0) + (n - 1.0) * a**(1.0 - 2.0 * n)


def simulate_bounce(w: float = 1.0 / 3.0, a_init: float = 3.0,
                    t_max: float = 12.0, n_points: int = 2000) -> BounceSolution:
    """Integrate a collapsing spin-fluid universe through the torsion bounce.

    Starts at a = a_init in the collapsing branch (H < 0), set by the Friedmann
    first integral, and integrates forward through the bounce into re-expansion.
    """
    n = _n(w)

    def H2(a):
        return a**(-n) * (1.0 - a**(-n))

    v_init = -a_init * math.sqrt(max(H2(a_init), 0.0))   # collapsing root

    def rhs(t, y):
        a, v = y
        return [v, _accel(a, n)]

    sol = solve_ivp(rhs, (0.0, t_max), [a_init, v_init],
                    t_eval=np.linspace(0.0, t_max, n_points),
                    rtol=1e-10, atol=1e-12, method="DOP853")
    a = sol.y[0]
    v = sol.y[1]
    H = v / a
    rho = a**(-n)
    i_min = int(np.argmin(a))
    return BounceSolution(t=sol.t - sol.t[i_min], a=a, H=H, rho=rho,
                          a_min=float(a[i_min]), rho_max=float(rho[i_min]), w=w)


def simulate_gr_collapse(w: float = 1.0 / 3.0, a_init: float = 3.0,
                         n_points: int = 2000) -> dict:
    """The same collapse in plain GR (no spin term): hits a = 0 in finite time.

    H^2 = (8 pi G/3) rho = a^{-n}, so da/dt = -a^{1 - n/2}; analytic until a -> 0.
    """
    n = _n(w)
    # da/dt = -a^{1-n/2}  ->  a(t) = (a_init^{n/2} - (n/2) t)^{2/n}
    t_sing = a_init**(n / 2.0) / (n / 2.0)
    t = np.linspace(0.0, t_sing * 0.999, n_points)
    a = (a_init**(n / 2.0) - (n / 2.0) * t)**(2.0 / n)
    return {"t": t, "a": a, "t_singularity": float(t_sing)}


def physical_timescale(rho_C: float = 1.0e50) -> float:
    """Dimensionless time unit in seconds: 1 / sqrt(8 pi G rho_C / 3)."""
    return 1.0 / math.sqrt(8.0 * math.pi * k.G * rho_C / 3.0)


def bounce_fwhm(sol: BounceSolution) -> float:
    """Proper-time full width of the density peak at half its maximum.

    Exact closed form: rho/rho_C = 1/2 at a = 2^{1/n}, and integrating
    dtau = da/(a|H|) from the bounce out to that radius gives, with n = 3(1+w),

        FWHM = (2/n) * integral_{1/2}^{1} u^{-3/2}(1-u)^{-1/2} du = 4/n.

    (The leading small-amplitude estimate 2 sqrt(2)/n undershoots by sqrt(2).)
    """
    rho = sol.rho
    t = sol.t
    i = int(np.argmin(sol.a))           # density peak
    half = sol.rho_max / 2.0

    def cross(idxs):
        # linear interpolation of the t where rho crosses `half`
        for a, b in zip(idxs[:-1], idxs[1:]):
            if (rho[a] - half) * (rho[b] - half) <= 0:
                f = (half - rho[a]) / (rho[b] - rho[a])
                return t[a] + f * (t[b] - t[a])
        return math.nan

    t_left = cross(range(i, -1, -1))
    t_right = cross(range(i, len(t)))
    return float(t_right - t_left)
