r"""The supraverse census: tabulate viable universes by type, and place our footprint.

The wallpaper math (population.py, accretion.py, genesis.py) lets us count the viable-
observer population along three independent axes and read off where we sit -- our
"footprint" -- plus the observable tags that go with each cell.

  GENERATION depth n (BHU_n): geometric, P(n | n>=1) = (1-eps) eps^{n-1} (population.py).
    Set by the DM+DE anchor (n>=1) + the mass-budget branching tail. NOT set by spin.
  CHIRALITY: matter vs antimatter, 50/50 by CPT (Ch.4). A lineage is pure; the two
    halves are observationally identical (a CPT mirror).
  FAMILY: clean/photon-dominated (eta ~ 1ppb, fed by a fair-sample/horizon-scale hole)
    vs baryon-rich/degenerate (eta -> 1, fed by concentrated baryons) (accretion.py).
    Observers like us are in the clean family.

Our footprint, with the tags we have derived/observe:
  * generation : n = 1 most likely (~90% at eps~0.1) -- we are a BHU1.
  * chirality  : matter (by convention; our half).
  * family     : clean (eta ~ 6e-10 ~ eta_cosmic -> fair-sample progenitor).
  * spin       : LOW/quiet is typical (viable seeds cluster just above threshold), but
    spin is LOCAL and re-drawn each generation -- it tags the progenitor hole, not n.
  * size       : ~a Hubble mass (Nariai-capped), ~M_vis.

IMPORTANT (corrects a tempting error): a shared-axis NULL does NOT promote us from
BHU1 to BHU2. Spin and generation are orthogonal -- spin is the local progenitor
hole's Kerr spin, generation is set by the anchor + branching. A quiet axis just means
our progenitor hole was low-spin; it leaves the generation census untouched.

The census is the framework's "where are we in the population" statement, made into a
table; matching our derived tags to the dominant cell is the footprint check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import population as pop


@dataclass(frozen=True)
class Footprint:
    generation: int
    chirality: str       # 'matter' | 'antimatter'
    family: str          # 'clean' | 'baryon-rich'
    fraction: float      # share of viable observers in this cell


def generation_pmf(epsilon: float = 0.1, n_max: int = 6):
    """P(BHU_n | n>=1) = (1-eps) eps^{n-1}: the geometric generation census."""
    n, p = pop.generational_depth_pmf(epsilon, n_max=n_max)
    return n, p


def chirality_split() -> dict:
    """50/50 matter/antimatter by CPT; the two are observationally identical."""
    return {"matter": 0.5, "antimatter": 0.5}


def family_split(f_clean: float = 0.9) -> dict:
    """Fraction of viable observers in the clean (photon-dominated) family vs the
    baryon-rich (degenerate) family. Viable, star-forming, observer-bearing universes
    overwhelmingly sit in the clean family (eta ~ ppb), the baryon-rich tail being
    degenerate; f_clean ~ 0.9 is the fiducial split."""
    return {"clean": f_clean, "baryon-rich": 1.0 - f_clean}


def census(epsilon: float = 0.1, f_clean: float = 0.9, n_max: int = 4):
    """Full (generation x chirality x family) table of viable-observer fractions.
    Returns a list of (generation, chirality, family, fraction), summing to 1."""
    n, pgen = generation_pmf(epsilon, n_max=n_max)
    chir = chirality_split()
    fam = family_split(f_clean)
    rows = []
    for ni, pn in zip(n, pgen):
        for cname, cp in chir.items():
            for fname, fp in fam.items():
                rows.append((int(ni), cname, fname, float(pn * cp * fp)))
    return rows


def our_footprint(epsilon: float = 0.1, f_clean: float = 0.9) -> Footprint:
    """The cell we occupy: BHU1, matter (our half), clean family -- and its share of
    all viable observers."""
    n, pgen = generation_pmf(epsilon)
    frac = pgen[0] * 0.5 * f_clean
    return Footprint(generation=1, chirality="matter", family="clean",
                     fraction=float(frac))


def footprint_is_typical(epsilon: float = 0.1, f_clean: float = 0.9,
                         threshold: float = 0.1) -> bool:
    """Is our cell a typical (dominant) one? True if it holds more than `threshold`
    of all viable observers -- i.e. we are where most observers are (Copernican)."""
    return our_footprint(epsilon, f_clean).fraction > threshold


def spin_changes_generation() -> bool:
    """A shared-axis result does NOT change our generation -- spin is local, n is set
    by the anchor + branching. (Guard against the spin->generation conflation.)"""
    return False
