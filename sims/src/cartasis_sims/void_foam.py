r"""The condensed void as a nucleation-and-growth foam: the average OGU size.

Individual OGUs grow steadily (rate-limited at the causal speed, never depleting
the infinite void on their own; growth.py). But at the POPULATION level there is a
shared-resource limit: OGUs nucleate throughout the void and grow until they
\emph{impinge} on neighbours -- the void between them is consumed. That collective
"running out of nothing" sets a characteristic cell size and hence an average OGU
mass, fixed by the one number we do not yet know: the birth (nucleation) rate.

This is exactly the Kolmogorov--Johnson--Mehl--Avrami (KJMA) nucleation-and-growth
problem. OGUs nucleate at rate beta (per comoving 3-volume per time) and grow a
horizon of radius R = c (t - t_birth). The extended (Avrami) volume fraction is

    X_ext(t) = beta * Integral_0^t (4/3) pi [c(t-t')]^3 dt' = (pi/3) beta c^3 t^4,

so the void fills (X_ext ~ 1) at

    t_fill = (3 / (pi beta c^3))^{1/4},

the characteristic cell radius is L ~ c t_fill ~ (c/beta)^{1/4}, and -- since a
filled cell is a horizon-scale hole of that radius -- the average OGU mass is

    M_avg ~ c^3 t_fill / (2 G).

Frequent births (large beta) tile the void into many small universes; rare births
let each grow large before impinging. So the average OGU size is the birth rate
read through the growth law -- the same single unknown (Q12) that the OGU mass and
the supraverse age all reduce to (genesis.py). A 2D Johnson--Mehl tessellation
(johnson_mehl_2d) visualises the foam and its cell-size spread; note the impinged
cells are curved polygons, not circles -- only an isolated, not-yet-impinged OGU
is round.
"""

from __future__ import annotations

import numpy as np

from . import constants as k


def fill_time(beta: float) -> float:
    """KJMA fill time t_fill = (3/(pi beta c^3))^{1/4} [s] for nucleation rate beta
    [per m^3 per s] and causal growth speed c."""
    return (3.0 / (np.pi * beta * k.c**3)) ** 0.25


def char_size(beta: float) -> float:
    """Characteristic cell (OGU) radius L ~ c t_fill [m]."""
    return k.c * fill_time(beta)


def avg_ogu_mass(beta: float) -> float:
    """Average OGU mass: horizon mass at impingement, c^3 t_fill/(2G) [kg]."""
    return k.c**3 * fill_time(beta) / (2.0 * k.G)


def number_density(beta: float) -> float:
    """OGUs per comoving 3-volume at impingement, ~ 1/L^3 [m^-3]."""
    return 1.0 / char_size(beta) ** 3


def birth_rate_for_mass(M_avg: float) -> float:
    """Invert avg_ogu_mass: the birth rate implied by an average OGU mass."""
    t = 2.0 * k.G * M_avg / k.c**3
    return 3.0 / (np.pi * k.c**3 * t**4)


def johnson_mehl_2d(n_seeds: int = 60, grid: int = 240, t_max: float = 1.0,
                    speed: float = 1.0, seed: int = 0):
    """2D Johnson--Mehl tessellation of the foam (unit box, periodic-free).

    Seeds nucleate at uniform random positions and birth times in [0, t_max]; a
    point is claimed by the seed minimising (t_birth + distance/speed) -- the first
    growth front to reach it. Returns (labels[grid,grid], areas[n_seeds]) with
    areas the fractional cell areas (phantom, overgrown seeds get ~0)."""
    rng = np.random.default_rng(seed)
    xs = rng.random(n_seeds)
    ys = rng.random(n_seeds)
    ts = rng.random(n_seeds) * t_max
    gx = (np.arange(grid) + 0.5) / grid
    X, Y = np.meshgrid(gx, gx)
    arrival = np.empty((n_seeds, grid, grid))
    for i in range(n_seeds):
        d = np.hypot(X - xs[i], Y - ys[i])
        arrival[i] = ts[i] + d / speed
    labels = np.argmin(arrival, axis=0)
    areas = np.array([(labels == i).mean() for i in range(n_seeds)])
    return labels, areas, (xs, ys, ts)
