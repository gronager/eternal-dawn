r"""Generations from the electroweak S budget -- the cap, and why it requires walking.

The idea (your conjecture): the number of generations is capped because each adds to
the electroweak S, and the total composite S must stay under the precision bound
(~0.1). So N_max ~ budget / S_per_generation.

This module makes the mechanism and its SENSITIVITY explicit. Two robust facts and one
sensitive number:

ROBUST 1 -- a cap EXISTS. If every generation contributes a roughly equal S, the
budget forces a finite number. The Standard Model does not explain why the number of
generations is finite; here it is a consequence.

ROBUST 2 -- our existence DEMANDS walking. The leading-order (constituent-loop)
S_per_gen = N_c/6pi ~ 0.16 gives N_max = floor(0.1/0.16) = 0 -- it would forbid even
the one generation we are made of. So for the framework to be consistent with us
existing at all, S_per_gen must be FAR below leading order: the walking is not optional
decoration, it is required by our own existence. And the Fierz (G_A=G_V) + RPA showed
the torsion interaction supplies walking, in the right direction.

SENSITIVE -- the exact cap. N_max = budget/S_per_gen is a steep 1/x law:
  S_per_gen ~ 0.033 -> cap 3;  0.05 -> 2;  0.025 -> 4;  0.016 -> 6.
So cap = 3 needs S_per_gen in a narrow window [0.025, 0.033], a ~5x walk below leading.
Whether the torsion walks to exactly there (cap 3) vs nearby (cap 2 or 4) is the
highly-sensitive, still-uncomputed number -- it needs the full chiral RPA (lattice).

CONCLUSION: the framework predicts a FINITE number of generations and REQUIRES walking
to allow our three; three is squarely in the plausible window; but the exact count is
not pinned, and saying it is would be dishonest. The cap is real; the number is owed.
"""

from __future__ import annotations

import numpy as np

S_BUDGET = 0.1                       # rough electroweak-precision ceiling on total S
S_LEADING = 3.0 / (6.0 * np.pi)      # constituent-loop S per generation ~ 0.16


def n_max(s_per_gen, budget=S_BUDGET):
    """The generation cap: how many generations fit under the S budget."""
    return int(budget / s_per_gen)


def s_per_gen_for_cap(cap, budget=S_BUDGET):
    """The S-per-generation window that yields a given cap: (lower, upper]."""
    return budget / (cap + 1), budget / cap


def requires_walking(budget=S_BUDGET):
    """Leading-order S forbids our existence -> walking is mandatory. Returns the
    walk-down factor needed merely to allow the 3 generations we observe."""
    s_needed = budget / 3.0
    return {"s_leading": S_LEADING, "s_needed_for_3": s_needed,
            "walk_factor_needed": S_LEADING / s_needed,
            "leading_cap": n_max(S_LEADING)}
