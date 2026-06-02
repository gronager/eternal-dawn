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


def ogu_evaporation_lifetime(M: float) -> float:
    """Hawking evaporation time of an OGU, 5120 pi G^2 M^3/(hbar c^4) [s]. For
    OGU masses (1e53-1e65 kg) this is ~1e143-1e178 s -- astronomically longer than
    any supraverse timescale, so OGUs do not die by evaporation; they only coarsen
    by merging or by the radiative siphoning below."""
    return 5120.0 * np.pi * k.G**2 * M**3 / (k.hbar * k.c**4)


# ---------------------------------------------------------------------------
# Radiative coarsening of the packed foam (Hawking "siphoning")
# ---------------------------------------------------------------------------
# Once the void is packed there is no fresh void to accrete, so each OGU sits in a
# bath set by its neighbours' Hawking radiation (and the void's own de Sitter
# horizon temperature -- both ~ T_H of a horizon-mass hole ~ 1e-30 K; Phase 0,
# Ch.2). At that bath each OGU is right at M ~ M_crit: marginal balance. But a black
# hole has NEGATIVE heat capacity (T_H ~ 1/M), so the balance is UNSTABLE: an OGU
# slightly bigger than its neighbours is slightly colder, absorbs more than it
# emits, and grows; a slightly smaller one is hotter, emits more, and shrinks --
# shrinking makes it hotter still, a runaway to evaporation. So the foam coarsens by
# RADIATIVE SIPHONING: big (cold) OGUs eat the radiation of small (hot) ones, which
# evaporate away into their larger neighbours. This is Ostwald ripening for
# universe-holes, and it runs on the Hawking timescale (~1e143 s) -- immortal on any
# ordinary clock, but not eternal.
#
# Mean field: with Hawking power P_i ~ M_i^-2 and the bath fixed by total-mass
# conservation, the equivalent bath mass is M_bar = (sum M_i^2 / sum M_i^-2)^{1/4}
# and dM_i/dtau = M_i^-2 [ (M_i/M_bar)^4 - 1 ] in units of the Hawking time of the
# mean mass. M_i > M_bar grow, M_i < M_bar evaporate; total mass is conserved.

def foam_bath_mass(masses):
    """Conservation-defined bath equivalent mass M_bar = (sum M^2/sum M^-2)^{1/4}.
    An OGU at M = M_bar is in (unstable) radiative balance with the foam."""
    m = np.asarray(masses, dtype=float)
    return (np.sum(m**2) / np.sum(m**-2)) ** 0.25


def coarsening_rate(masses, m_dead: float = 0.02):
    """dM_i/dtau = M_i^-2 [(M_i/M_bar)^4 - 1] (Hawking-time units, mean=1); the
    negative-heat-capacity instability that drives radiative siphoning. Fully
    evaporated OGUs (M < m_dead) drop out of the bath and freeze, so the surviving
    foam's mass is conserved and the integrator avoids the M->0 stiffness."""
    m = np.asarray(masses, dtype=float)
    alive = m > m_dead
    rate = np.zeros_like(m)
    if alive.sum() == 0:
        return rate
    ma = m[alive]
    Mb = (np.sum(ma**2) / np.sum(ma**-2)) ** 0.25
    rate[alive] = ma**-2 * ((ma / Mb) ** 4 - 1.0)
    return rate


def simulate_coarsening(masses0, tau_max: float = 4.0, n_points: int = 400):
    """Evolve a foam of OGUs by radiative siphoning. masses0 in units of their
    initial mean; tau in Hawking times of the mean. Returns (tau, M[i, t])."""
    from scipy.integrate import solve_ivp

    def rhs(t, y):
        return coarsening_rate(y)

    sol = solve_ivp(rhs, (0.0, tau_max), np.asarray(masses0, dtype=float),
                    t_eval=np.linspace(0.0, tau_max, n_points),
                    rtol=1e-7, atol=1e-9, max_step=0.01)
    return sol.t, sol.y


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
