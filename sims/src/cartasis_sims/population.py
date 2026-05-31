r"""Supraverse population dynamics: birth vs growth, and our generation depth.

This module does, properly, what Chapter 5 previously only postulated: derive the
*distribution over generation depth* of universes in the supraverse foam, and ask
where our universe sits in it. The structure is a bounded branching process.

The model (minimal, stated assumptions)
---------------------------------------
1. Birth. Original-generation universes (OGUs == BHU_0) nucleate directly from the
   void at a steady rate per unit four-volume (the void is infinite and
   stationary). Their number born per unit time is a constant we can set to 1
   without loss of generality -- only *ratios* between generations matter for the
   depth distribution.

2. Growth / branching. Once a universe is viable it grows by runaway accretion
   (negative heat capacity above M_crit) and internally forms black holes. Each
   internal hole that (a) exceeds its internal-bath M_crit and (b) clears the
   chiral-vortical birth filter (|omega/T| > eta_min/C) becomes a viable child
   one generation deeper. Let N = mean number of *viable* children per universe
   (the effective branching ratio). N is NOT the ~1e18 holes per universe; it is
   the far smaller count that are both super-critical and pass the birth filter.

3. Truncation. The tree cannot be infinitely deep: each generation must retain
   enough structure-forming capacity to make its own holes. Let p = per-generation
   probability that a child remains viable (fertility-preserving). Equivalently a
   universe has, on average, m = N*p viable grandchildren-capable lines. We model
   truncation either as a hard depth cap D or, more honestly, as a sub-critical
   branching tail with effective reproduction m <= 1 that dies out.

Two regimes, two very different answers for "where are we":
  * m > 1 (supercritical): generation count grows like m^n; counting *universes*,
    the typical universe is as deep as truncation allows -> we are deep, n large.
  * m < 1 (subcritical): the lineage dies out; most universes are shallow.
  * m = 1 (critical): scale-free; broad distribution of depths.

The robust observational anchor
--------------------------------
Independently of the branching parameters, our universe HAS dark matter and dark
energy. In this framework both require a parent: dark matter is parent material
projected through the membrane, dark energy is the parent's ongoing accretion
(Ch. 6). An OGU (n=0), born from the void with no parent, would have neither.

    Therefore we are NOT an OGU. Our generation n >= 1.  (robust)

This is the one thing the data settle cleanly. The branching model then tells us
whether, given n>=1, we are likely shallow (n=1,2) or deep (n>>1).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Generation-depth distribution of a bounded / branching population
# ---------------------------------------------------------------------------
def generation_counts(m, D):
    """Relative number of universes at each generation n=0..D for branching
    ratio m (mean viable children per universe), counting universes (not
    observers). N_n proportional to m^n up to the truncation depth D.
    """
    n = np.arange(D + 1)
    return m ** n


def generation_pmf(m, D):
    """Normalized P(n): fraction of all universes at generation n (0..D)."""
    c = generation_counts(m, D)
    return c / c.sum()


def expected_generation(m, D):
    """Mean generation depth <n> over the universe-counted population."""
    n = np.arange(D + 1)
    p = generation_pmf(m, D)
    return float((n * p).sum())


def observer_weighted_pmf(m, D, structure_decay=1.0):
    r"""P(n) weighted by observers, not universes.

    If each generation passes a fraction `structure_decay` of the
    structure-forming capacity (galaxies/observers per universe) to the next,
    observer weight per generation-n universe ~ structure_decay^n. Combined with
    m^n universes, the observer share at depth n ~ (m * structure_decay)^n.
    structure_decay=1 reproduces the universe-counted pmf; <1 pulls the typical
    OBSERVER shallower than the typical universe.
    """
    n = np.arange(D + 1)
    w = (m * structure_decay) ** n
    return w / w.sum()


def survival_truncation_depth(p_viable, eps=1e-3):
    """Effective max depth where a fertility-preserving fraction p_viable per
    generation leaves >= eps of lineages still viable: p_viable^D = eps."""
    if not (0 < p_viable < 1):
        return np.inf
    return int(np.ceil(np.log(eps) / np.log(p_viable)))


# ---------------------------------------------------------------------------
# The robust observational anchor: DM + DE => we have a parent => n >= 1
# ---------------------------------------------------------------------------
def we_are_ogu(has_dark_matter, has_dark_energy):
    """An OGU (n=0) has no parent, hence no projected dark matter and no
    injection-driven dark energy. Returns True only if BOTH are absent."""
    return (not has_dark_matter) and (not has_dark_energy)


def min_generation_from_observations(has_dark_matter=True, has_dark_energy=True):
    """Lower bound on our generation from what we observe. Presence of either
    parent-sourced component forces n >= 1."""
    return 1 if (has_dark_matter or has_dark_energy) else 0


# ---------------------------------------------------------------------------
# Posterior on our depth: branching prior x observational floor
# ---------------------------------------------------------------------------
def depth_posterior(m, D, n_min=1, structure_decay=1.0):
    """P(n | branching m,D AND n>=n_min), observer-weighted. Renormalize the
    observer-weighted pmf after zeroing generations below the observational
    floor (we know we are not an OGU)."""
    p = observer_weighted_pmf(m, D, structure_decay)
    p = p.copy()
    p[:n_min] = 0.0
    s = p.sum()
    return p / s if s > 0 else p


# ---------------------------------------------------------------------------
# Explicit read-off: are we likely BHU1 or BHU2?
# ---------------------------------------------------------------------------
# Counting universes, generation n holds N_n ~ m^n of them, so given the robust
# floor n >= 1 (we are not an OGU), P(BHU_n) ~ m^n / sum. For a SUBCRITICAL
# branching m < 1 this is geometric, P(BHU_n) = (1-m) m^{n-1}: lineages die out,
# the population is shallow, and we are most likely BHU1. For SUPERCRITICAL m > 1
# the population is dominated by the deepest generation the supraverse has had
# TIME to grow (truncation D), so we are almost surely deep and being BHU1-2 is
# exponentially unlikely. The "are we shallow?" question therefore reduces to
# whether the effective viable-children-per-universe m is below 1 -- or,
# equivalently, whether the supraverse is too young (small D) to have grown deep.

def prob_bhu(m, n, D=400):
    """P(we are BHU_n | n >= 1) for branching ratio m and truncation depth D."""
    D = min(int(D), 100000)              # guard against huge-D array allocation
    ns = np.arange(1, D + 1, dtype=float)
    w = m ** ns
    w = w / w.sum()
    return float(w[n - 1]) if 1 <= n <= D else 0.0


def shallow_probability(m, D=400, n_max=2):
    """P(we are within the first n_max generations | n >= 1): the chance we are
    'BHU1 or BHU2' (n_max=2). Near 1 for strongly subcritical m, ~0 for m > 1."""
    return float(sum(prob_bhu(m, n, D) for n in range(1, n_max + 1)))


# ---------------------------------------------------------------------------
# OGU size (mass) distribution from birth vs growth
# ---------------------------------------------------------------------------
# Births inject OGUs at small mass: the Cartasis-density fluctuation probability
# falls as exp(-M c^2 tau/hbar), so dN_birth/dM ~ exp(-M/M0) with M0 = hbar/(c^2
# tau) a small birth scale. Once born, super-critical holes grow by runaway
# accretion at rate Mdot(M) = A M^g (g=1 Eddington-like, g=2 Bondi-like). In a
# stationary foam the number density n(M) obeys a continuity equation in mass
# space, d/dM[Mdot(M) n(M)] = B(M) - death. For M well above the birth scale the
# upward flux Mdot*n equals the (nearly constant) total birth rate, so
#
#     n(M) ~ 1 / Mdot(M) ~ M^{-g}              (between M0 and a death cutoff),
#
# i.e. a POWER LAW: many small universes, few large ones, derived rather than
# assumed. A high-mass cutoff M_max comes from environment depletion / the
# stationary balance; a low-mass VIABILITY cut M_vis (enough mass for astrophysics
# and black-hole formation) selects the universes that can host observers and seed
# descendants. Because n(M) falls with M, observer-bearing universes pile up just
# above M_vis -- we should sit near the viability edge, which is the quantitative
# content of the old "we are near the optimum" conjecture.

def ogu_birth_rate(M, M0=1.0):
    """Birth rate per unit mass: B(M) ~ exp(-M/M0), small M favoured."""
    M = np.asarray(M, dtype=float)
    return np.exp(-M / M0)


def ogu_mass_density(M, g=2.0, M0=1.0, M_max=1e3):
    r"""Stationary OGU number density n(M) ~ M^{-g} exp(-M/M_max), valid above the
    birth scale M0. g is the growth-law exponent (Mdot ~ M^g): g=1 Eddington,
    g=2 Bondi. Returns an unnormalized density (use on a log-M grid)."""
    M = np.asarray(M, dtype=float)
    return M ** (-g) * np.exp(-M / M_max) * (1.0 - np.exp(-M / M0))


def observer_mass_density(M, g=2.0, M0=1.0, M_max=1e3, M_vis=10.0):
    """OGU mass density restricted to viable (observer-bearing) universes:
    zero below the viability mass M_vis, n(M)~M^{-g} above. The peak of this is
    the most probable host-universe mass -> our expected neighbourhood."""
    M = np.asarray(M, dtype=float)
    dens = ogu_mass_density(M, g, M0, M_max)
    return np.where(M >= M_vis, dens, 0.0)


# ---------------------------------------------------------------------------
# Spin (vorticity) distribution: birth favours low spin, viability needs spin
# ---------------------------------------------------------------------------
# The OG seed's net vorticity omega sets its baryon asymmetry, eta ~ C|omega|/T
# (Ch. 4): the SIGN of omega fixes matter vs antimatter (CPT-symmetric, so the
# birth distribution is even in omega -> half matter, half antimatter), and the
# MAGNITUDE fixes purity. A random fluctuation most likely has LOW net vorticity,
# P_birth(omega) ~ Gaussian(0, sigma) peaked at 0. But a universe with
# |eta| < eta_min cannot complete baryogenesis: matter and antimatter annihilate
# to near-symmetric radiation, no structure, no black holes -- a "hellish",
# sterile universe that simply evaporates and seeds nothing. Viability therefore
# requires |omega| > omega_min = eta_min T / C. Among viable seeds, higher spin
# means higher purity and more descendants (productivity), but is exponentially
# rarer. The observer/seed-weighted spin distribution is thus
#
#     P_obs(omega) ~ P_birth(omega) * [productivity(omega)] * 1[|omega|>omega_min],
#
# which peaks just ABOVE the viability threshold: viable universes are typically
# LOW-purity, just past sterile. Our small observed asymmetry (eta ~ 1e-9..1e-8)
# places us near that edge -- a prediction testable against any net rotation of
# our universe (galaxy-spin handedness, CMB axis hints).

def spin_birth_pdf(omega, sigma=1.0):
    """Birth vorticity distribution: Gaussian peaked at 0 (low spin likeliest)."""
    omega = np.asarray(omega, dtype=float)
    return np.exp(-0.5 * (omega / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def spin_threshold(eta_min, C, T):
    """Minimum viable vorticity for baryogenesis: omega_min = eta_min * T / C."""
    return eta_min * T / C


def purity_from_spin(omega, C, T):
    """Baryon asymmetry (purity) from vorticity: eta = C|omega|/T."""
    return C * np.abs(np.asarray(omega, dtype=float)) / T


def spin_observer_pdf(omega, sigma=1.0, omega_min=0.5, prod_power=1.0):
    r"""Seed-weighted spin distribution: P_birth * productivity * viability gate.
    productivity ~ (|omega|-omega_min)^prod_power above threshold (more spin ->
    more purity -> more structure and descendants), zero below. Returns an
    unnormalized density; normalize on a grid for plotting."""
    omega = np.asarray(omega, dtype=float)
    prod = np.where(np.abs(omega) > omega_min,
                    (np.abs(omega) - omega_min) ** prod_power, 0.0)
    return spin_birth_pdf(omega, sigma) * prod


def sterile_fraction(sigma=1.0, omega_min=0.5):
    """Fraction of OG seeds that are sub-threshold (sterile, evaporate): the
    Gaussian mass within |omega| < omega_min."""
    from math import erf
    return float(erf(omega_min / (sigma * np.sqrt(2.0))))

