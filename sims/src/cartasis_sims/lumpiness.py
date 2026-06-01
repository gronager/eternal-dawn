r"""The lumpiness ladder: inhomogeneity accumulating across Cartasis membranes.

A structural prediction of the recursive supraverse, and a quantitative reading of one
of LambdaCDM's real tensions -- the excess of ULTRA-LARGE structures (the Giant Arc,
the Big Ring, the Sloan Great Wall, the KBC void), which sit uncomfortably with the
cosmological principle.

The picture:
  * The VOID is homogeneous (lumpiness L = 0) but unstable to quantum fluctuations;
    rare ones nucleate OGUs -- the first structure.
  * An OGU has NO parent, so NO inherited lumpiness; it makes only its own (weak,
    dark-matter-free) internal structure -- minimal large-scale power. THIS is what
    LambdaCDM effectively assumes: a smooth/singular start, building structure only
    internally. In SCT terms, LambdaCDM predicts an OGU.
  * A BHU_n INHERITS its parent's lumpiness (the lumpy meal projected through the
    membrane, Ch.6) AND adds its own internal structure on top. So large-scale power
    ACCUMULATES generation by generation across membranes.

A minimal model: write the dimensionless large-scale lumpiness L_n (excess power on
scales above a universe's own homogeneity scale) as a quadrature sum of an intrinsic
per-generation contribution a and a transferred fraction b of the parent's:

    L_n^2 = a^2 + (b L_{n-1})^2,   L_0 (OGU) = small,

which converges to a fixed point L* = a / sqrt(1 - b^2) for b < 1. So:
  L(void)=0 < L(OGU) small < L(BHU1) < L(BHU2) < ... -> L*,
a monotonic LADDER of increasing large-scale lumpiness as one dives through membranes.
The jump from OGU to BHU1 is the largest (the parent's lumpiness first appears); deeper
generations add diminishing increments toward L*.

The observable: our universe shows MORE ultra-large structure than a smooth-start
(OGU-like, LambdaCDM) model predicts -- the Giant Arc (~1 Gly) and Big Ring (~0.4 Gpc)
exceed the ~260 Mpc homogeneity scale and are hard to assemble internally in the
available time. In SCT this excess is the inherited rung of the ladder: evidence that
we are a BHU (n>=1), not an OGU. It is yet another generation-depth indicator,
alongside dark matter and dark energy. (The model is illustrative -- a and b are not
yet derived from membrane projection -- but the DIRECTION, monotonic excess large-scale
power for n>=1, is the robust, falsifiable content.)
"""

from __future__ import annotations

import numpy as np

MPC_PER_GLY = 306.6
HOMOGENEITY_MPC = 260.0          # ~ the LambdaCDM homogeneity scale (End of Greatness)


def lumpiness_ladder(a: float = 1.0, b: float = 0.7, L0: float = 0.3,
                     n_max: int = 6) -> np.ndarray:
    """L_n for n=0..n_max via L_n^2 = a^2 + (b L_{n-1})^2, L_0 = OGU lumpiness.
    Returns the array [L_OGU, L_BHU1, L_BHU2, ...]."""
    L = [float(L0)]
    for _ in range(n_max):
        L.append(float(np.sqrt(a**2 + (b * L[-1]) ** 2)))
    return np.array(L)


def fixed_point(a: float = 1.0, b: float = 0.7) -> float:
    """Asymptotic lumpiness L* = a/sqrt(1-b^2) for b<1 (deep generations saturate)."""
    if not (0.0 <= b < 1.0):
        return float("inf")
    return a / np.sqrt(1.0 - b**2)


def is_monotonic(L: np.ndarray) -> bool:
    """The ladder increases with generation depth (more lumpiness deeper in)."""
    return bool(np.all(np.diff(L) > 0))


def biggest_jump_index(L: np.ndarray) -> int:
    """Which membrane crossing adds the most lumpiness (the OGU->BHU1 step)."""
    return int(np.argmax(np.diff(L)))


def gly_to_k(gly: float) -> float:
    """Comoving scale in Gly -> wavenumber k [1/Mpc] (k = 2 pi / lambda)."""
    return 2.0 * np.pi / (gly * MPC_PER_GLY)


def exceeds_homogeneity(scale_mpc: float) -> bool:
    """True if a structure is larger than the homogeneity scale -- i.e. it is
    'too big' for a smooth-start (OGU/LambdaCDM) model to assemble comfortably."""
    return scale_mpc > HOMOGENEITY_MPC


# representative ultra-large structures (comoving extent)
ULTRA_LARGE = {
    "BAO scale (ref)":      150.0,
    "Sloan Great Wall":     420.0,
    "Giant Arc":            1000.0 * 306.6 / 1000.0 * 1.0,   # ~1 Gly extent -> Mpc
    "Big Ring (diameter)":  400.0,
    "KBC void (radius)":    300.0,
}
# fix Giant Arc to its ~1 Gly extent in Mpc
ULTRA_LARGE["Giant Arc"] = 1.0 * MPC_PER_GLY
