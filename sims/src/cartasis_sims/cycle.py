r"""The recursive cycle: a heat-dead universe is the void of the next foam.

The far future of any universe is dark-energy-dominated de Sitter: everything
redshifts away, black holes evaporate, and the interior becomes empty, cold
(T_dS ~ 1e-30 K), and -- crucially -- conformally SMOOTH. By Penrose's Weyl
Curvature Hypothesis the gravitational (clumping/Weyl) entropy of such a state is
LOW, even though the total entropy is maximal. So a heat-dead interior is exactly a
low-gravitational-entropy, CPT-symmetric, fluctuating quantum vacuum: the "void"
SCT nucleates OGUs from. The void is therefore not a primordial substrate outside
everything; it is the asymptotic state every universe reaches. A ginormous heat-dead
OGU breeds the next foam of OGUs inside itself, and it cycles -- forever, with no
first cause.

This is Penrose's CCC (the smooth far future = a fresh aeon's low-Weyl beginning),
made to BRANCH: a dead aeon does not hand off to one successor by conformal
rescaling but nucleates MANY OGUs in its emptiness, each seeding a lineage. SCT is
"CCC + foam".

The apparent paradox -- heat death is maximal entropy, yet nucleation needs low
entropy -- is resolved by Penrose's distinction. TOTAL entropy (thermal + horizon)
rises monotonically (the Second Law). The WEYL/clumping entropy cycles: low at the
bounce (smooth), high while structure and black holes form, low again once the holes
evaporate and the interior smooths to de Sitter. It is the Weyl entropy that a fresh
bounce needs low, and it resets every cycle. The reset runs on the black-hole
Hawking timescale (~1e100 yr) -- the same clock as the foam coarsening
(void_foam.py).

This module gives a toy model of the two entropies over a universe's life and the
de Sitter nucleation framing that ties the cycle's birth rate to dark energy.
"""

from __future__ import annotations

import numpy as np

from . import constants as k

# characteristic times in a universe's life [s]
T_BOUNCE = 1.0e-36           # ~ bounce / reheating
T_STRUCTURE = 1.0e16         # stars, galaxies, black holes form (~Gyr)
T_BH_EVAP = 1.0e108          # last black holes Hawking-evaporate (~1e100 yr)
T_DESITTER = 1.0e114         # interior smoothed to empty de Sitter (heat death)

# characteristic entropies [k_B units], order-of-magnitude (our universe)
LOG_S_RAD = 88.0             # radiation entropy of the observable universe
LOG_S_BH = 103.0             # entropy in black holes at peak structure
LOG_S_DESITTER = 122.0       # de Sitter cosmological-horizon entropy (the maximum)


def _up(t, tc, w=0.7):
    """Smooth log-time step from 0 to 1 around tc (width w decades)."""
    return 1.0 / (1.0 + 10.0 ** (-(np.log10(t) - np.log10(tc)) / w))


def weyl_entropy(t):
    """Gravitational (Weyl/clumping) entropy [log10 k_B]: low -> high (black holes)
    -> low (smooth de Sitter). This is what a fresh bounce needs small, and it
    RESETS each cycle. Modelled as a bump rising at structure formation and falling
    as the holes evaporate."""
    t = np.asarray(t, dtype=float)
    base = LOG_S_RAD - 8.0                       # nearly smooth floor
    bump = (LOG_S_BH - base) * _up(t, T_STRUCTURE, 0.8) * (1.0 - _up(t, T_BH_EVAP, 1.6))
    return base + bump


def total_entropy(t):
    """Total entropy [log10 k_B]: thermal + horizon, monotonically increasing to the
    de Sitter maximum (the Second Law). Black-hole entropy is a step on the way; it
    does not decrease when holes evaporate (it becomes radiation, then horizon)."""
    t = np.asarray(t, dtype=float)
    s = 10.0 ** LOG_S_RAD
    s = s + 10.0 ** LOG_S_BH * _up(t, T_STRUCTURE, 0.8)
    s = s + 10.0 ** LOG_S_DESITTER * _up(t, T_DESITTER, 2.0)
    return np.log10(s)


def de_sitter_temperature(H_lambda: float) -> float:
    """Gibbons-Hawking temperature of the de Sitter void, T_dS = hbar H/(2 pi kB) [K].
    For our Lambda (H ~ 1.8e-18 s^-1) this is ~1e-30 K -- the same ~T_H of a
    horizon-mass hole that sets the foam-coarsening bath (void_foam.py)."""
    return k.hbar * H_lambda / (2.0 * np.pi * k.kB)


def cycle_closes(nucleation_rate_per_volume_time: float = 0.0) -> bool:
    """The de Sitter void lasts forever, so even an exponentially tiny nucleation
    rate eventually fires: the cycle always closes, given infinite time."""
    return nucleation_rate_per_volume_time >= 0.0
