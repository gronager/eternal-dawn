r"""Growth of a universe-hole in the void, and the depth it implies.

Closing the population loop (population.py): birth gives a seed mass, growth
turns it into a universe and sets how fast generations stack, hence our likely
depth. The dynamics are nearly parameter-free once written through the Hawking
formula.

A hole of mass M in a bath at temperature T_bath emits Hawking radiation at
T_H = hbar c^3/(8 pi G M kB) and absorbs the bath with the same area/greybody
factor, so the NET power is the Hawking power times the temperature contrast:

    dM/dt = (P_H/c^2) [ (T_bath/T_H)^4 - 1 ],   P_H = hbar c^6/(15360 pi G^2 M^2).

Since T_H ~ 1/M, T_bath/T_H = M/M_crit with the critical mass

    M_crit = hbar c^3 / (8 pi G kB T_bath),

so in units x = M/M_crit and t' = t/tau0, tau0 = 3 t_Hawking(M_crit), the law is
universal:

    dx/dt' = (x^4 - 1) / x^2.

Below x = 1 (M < M_crit) the hole evaporates; above it the bath wins and the hole
runs away as dM/dt ~ M^2 (Bondi-like, g = 2 -- this is the exponent behind the
n(M) ~ M^{-2} OGU mass power law in population.py). The runaway LOOKS like a
finite-time blow-up, but it cannot blow up: CAUSALITY caps it. A hole cannot
accrete matter outside its causal horizon, so dM/dt saturates at the rate the
horizon mass itself grows,

    dM/dt <= c^3/(2G) ~ 2e35 kg/s ~ 3e12 Msun/yr  (causal_growth_rate),

and thereafter the hole grows STEADILY, M ~ c^3 t/(2G) (the horizon-mass law of
genesis.py). Whether it ever stops then depends on the environment, and the two
regulators part here:

  * Infinite void (an OGU): fresh vacuum perpetually crosses into the growing
    horizon, so it never runs out -- rate-limited and unbounded, M ~ c^3 t/(2G).
  * Finite parent (an internal hole): the parent's dark-energy expansion carries
    matter beyond reach faster than the hole can eat, so it DEPLETES its accessible
    patch and freezes at M_final ~ patch mass -- "runs out of nothing."

The depletion model below is the finite-environment (descendant) case; the causal
rate and genesis.py give the infinite-void (OGU) case. tau_gen ~ tau0 still sets
the generation time and the high-mass cutoff M_max for the finite case.

Depth. Each generation needs ~tau_gen to grow up and spawn the next, so a
supraverse of age T_supra has grown at most D_max = T_supra/tau_gen generations.
This is a HARD cap on our depth independent of the branching ratio m: if growth
is slow or the supraverse young (D_max <= 2) we are forced to be BHU1-2; if growth
is fast and the supraverse old (D_max large) the depth is set by m
(population.prob_bhu). The void bath temperature T_bath sets tau0 and hence
everything -- the one big unknown (Q15a, Q15c).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

from . import constants as k


def causal_growth_rate() -> float:
    """Max accretion rate set by causality: c^3/(2G) [kg/s]. A hole cannot eat
    matter outside its horizon, so dM/dt saturates here and M ~ c^3 t/(2G)."""
    return k.c**3 / (2.0 * k.G)


def m_crit(T_bath: float) -> float:
    """Critical mass: above this T_H < T_bath and the hole grows. [kg]"""
    return k.hbar * k.c**3 / (8.0 * math.pi * k.G * k.kB * T_bath)


def hawking_time(M: float) -> float:
    """Hawking evaporation time 5120 pi G^2 M^3/(hbar c^4). [s]"""
    return 5120.0 * math.pi * k.G**2 * M**3 / (k.hbar * k.c**4)


def tau0(T_bath: float) -> float:
    """Growth time unit: 3 x the Hawking time of M_crit. [s]"""
    return 3.0 * hawking_time(m_crit(T_bath))


def growth_rate_universal(x):
    """dx/dt' = (x^4 - 1)/x^2 in units x = M/M_crit, t' = t/tau0 (static bath)."""
    x = np.asarray(x, dtype=float)
    return (x**4 - 1.0) / x**2


def runaway_time_universal(x0: float) -> float:
    r"""Time t' for a static-bath hole to run away from x0 (x -> infinity).

    t' = \int_{x0}^\infty x^2/(x^4-1) dx, finite for x0 > 1 (~ 1/x0 for x0 >> 1).
    """
    if x0 <= 1.0:
        return math.inf
    from scipy.integrate import quad
    val, _ = quad(lambda x: x**2 / (x**4 - 1.0), x0, np.inf, limit=200)
    return float(val)


def simulate_runaway(x0: float = 2.0, t_max: float = 1.0, n: int = 2000):
    """Static-bath trajectory x(t'); blows up near runaway_time_universal(x0)."""
    def rhs(t, y):
        return [growth_rate_universal(y[0])]

    def blow(t, y):
        return y[0] - 1e6
    blow.terminal = True
    blow.direction = 1.0
    sol = solve_ivp(rhs, (0.0, t_max), [x0], events=blow, rtol=1e-9, atol=1e-12,
                    dense_output=True, max_step=1e-3)
    t = np.linspace(0.0, sol.t[-1], n)
    return t, sol.sol(t)[0]


def simulate_depletion(x0: float = 2.0, bath: float = 50.0,
                       t_max: float = 40.0, n: int = 4000):
    """Deplete a finite local bath: the hole eats it and freezes at ~x0+bath.

    db/dt' = -dx/dt' (mass conservation); once the bath is gone only Hawking
    remains. Returns (t', x, b)."""
    def rhs(t, y):
        x, b = y
        x = max(x, 1e-9)
        frac = max(b, 0.0) / bath
        dx = (x**4 * frac - 1.0) / x**2
        if b <= 0.0:
            dx = min(dx, 0.0)
        return [dx, -dx]

    sol = solve_ivp(rhs, (0.0, t_max), [x0, bath],
                    t_eval=np.linspace(0.0, t_max, n),
                    rtol=1e-8, atol=1e-11, max_step=1e-2)
    return sol.t, sol.y[0], sol.y[1]


def freeze_out_mass(x0: float = 2.0, bath: float = 50.0) -> float:
    """Final mass (units M_crit): the hole eats the patch, M_final ~ x0 + bath."""
    return x0 + bath


def generation_time(T_bath: float, x0: float = 2.0) -> float:
    """Physical time to grow up (run away) from a seed at x0 = M_seed/M_crit. [s]"""
    return runaway_time_universal(x0) * tau0(T_bath)


def max_depth(supra_age_s: float, T_bath: float, x0: float = 2.0) -> int:
    """Hard cap on generation depth: D_max = T_supra / tau_gen (>= 0)."""
    tg = generation_time(T_bath, x0)
    return int(supra_age_s // tg) if math.isfinite(tg) and tg > 0 else 0
