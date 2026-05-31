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
