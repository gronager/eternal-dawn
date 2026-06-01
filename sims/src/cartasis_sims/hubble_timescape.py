r"""The Hubble tension from inhomogeneity: the timescape clock-rate effect.

The H0 tension is a ~5-sigma gap between the early-universe value (CMB sound
horizon, H0 ~ 67.4) and the local distance-ladder value (Cepheid+SNe, H0 ~ 73.0).
SCT is natively an inhomogeneous, backreaction cosmology -- its dark sector is the
parent's ongoing accretion and its large-scale structure is a void-dominated foam
(Ch6) -- so it inherits Wiltshire's *timescape* resolution of the tension without
adding a parameter. This module makes that quantitative.

The mechanism (Wiltshire 2007, 2009). The real universe is a two-phase mix: VOIDS
(by volume fraction f_v, curvature-dominated, a_v ~ t) and WALLS (f_w = 1 - f_v,
where galaxies and observers sit, ~Einstein-de Sitter, a_w ~ t^{2/3}). General
relativity makes the clocks differ: a void clock runs FASTER than a wall clock
(less bound mass, larger volume per unit proper time). We are wall observers. The
quantity a wall observer infers -- the "dressed" Hubble parameter -- therefore
exceeds the volume-average "bare" Hubble parameter once voids dominate the volume.

Tracker solution. The exact Buchert two-phase tracker has the void fraction

    f_v(t) = t / (t + b),    b = (1 - f_v0)(2 + f_v0) / (3 f_v0 Hbar0),

with t0 = (2 + f_v0)/(3 Hbar0). So f_v -> 0 in the early universe (homogeneous:
the CMB epoch is UNAFFECTED) and grows to its observed f_v0 ~ 0.76 today.

Dressed vs bare Hubble (Wiltshire's closed form):

    H0_dressed = R(f_v0) Hbar0,   R(f_v) = (4 f_v^2 + f_v + 4) / (2 (2 + f_v)).

R(0) = 1 (no inhomogeneity, no effect); R rises through 1 as f_v grows, reaching
R(0.762) = 1.28 -- a 28% wall-vs-average Hubble excess. The present mean lapse
(clock-rate) is gamma-bar0 = 1.381 in the best fit; the dressed Hubble is slightly
below the raw lapse because it carries the -d(gamma-bar)/dt term.

The point for the tension. The effect has the RIGHT SIGN (local/wall > global) and
MORE than enough MAGNITUDE: at the observed void fraction it is 28%, against an
8.4% early-vs-late gap, and it switches on only at late times (f_v ~ 1e-4 at
recombination), so the CMB inference is untouched. The full data fit (Wiltshire et
al.) returns a single self-consistent dressed H0 ~ 61.7; what THIS module shows,
parameter-free, is that inhomogeneity at the *observed* void fraction is amply
capable of the tension in the right direction -- the "degrees of freedom" Ch8
claimed, made numerical. It is a sign-and-magnitude (scaling) result, not a
replacement for the global fit.

Reproducible companion to chapters/08-observational-tests.tex, "Hubble tension".
"""

from __future__ import annotations

# Anchors (km/s/Mpc).
H0_CMB = 67.4        # Planck 2018 early-universe (sound horizon)
H0_LOCAL = 73.04     # SH0ES 2022 local distance ladder (Riess et al.)

# Wiltshire timescape best-fit anchors.
FV0_BESTFIT = 0.762       # present void volume fraction
GAMMA0_BESTFIT = 1.381    # present mean phenomenological lapse (clock-rate)
H0_DRESSED_BESTFIT = 61.7  # dressed (wall-observer) Hubble, best fit


def dressed_over_bare(f_v: float) -> float:
    """Wiltshire's closed-form ratio of the dressed (wall-observer) Hubble to the
    bare (volume-average) Hubble, R(f_v) = (4 f_v^2 + f_v + 4)/(2(2 + f_v)).
    R(0) = 1 (homogeneous limit); R(0.762) ~ 1.28."""
    return (4.0 * f_v**2 + f_v + 4.0) / (2.0 * (2.0 + f_v))


def tracker_b(f_v0: float = FV0_BESTFIT, Hbar0: float = 1.0) -> float:
    """Tracker timescale b in f_v(t) = t/(t+b), in units of 1/Hbar0 when
    Hbar0=1. b = (1 - f_v0)(2 + f_v0)/(3 f_v0 Hbar0)."""
    return (1.0 - f_v0) * (2.0 + f_v0) / (3.0 * f_v0 * Hbar0)


def void_fraction(t_over_t0: float, f_v0: float = FV0_BESTFIT) -> float:
    """Tracker void fraction at time t (in units of t0). f_v(t) = t/(t+b) with
    t0 = (2+f_v0)/(3 Hbar0), so b/t0 = (1-f_v0)/f_v0. Returns f_v0 at t=t0,
    -> 0 as t -> 0 (homogeneous early universe)."""
    b_over_t0 = (1.0 - f_v0) / f_v0
    return t_over_t0 / (t_over_t0 + b_over_t0)


def void_fraction_at_recombination(z_rec: float = 1100.0,
                                   f_v0: float = FV0_BESTFIT) -> float:
    """Void fraction at recombination. In the matter era t ~ (1+z)^{-3/2}, so
    t/t0 ~ (1+z)^{-3/2}. Comes out ~1e-4: the clock effect is OFF at the CMB,
    which is why the early-universe H0 inference is unmodified."""
    t_over_t0 = (1.0 + z_rec) ** (-1.5)
    return void_fraction(t_over_t0, f_v0)


def wall_vs_average_offset(f_v0: float = FV0_BESTFIT) -> float:
    """Fractional excess of the dressed (wall/local) Hubble over the bare
    (volume-average/global) Hubble at the observed void fraction: R(f_v0) - 1.
    ~0.28 at f_v0 = 0.762."""
    return dressed_over_bare(f_v0) - 1.0


def observed_tension_fraction() -> float:
    """The measured early-vs-late gap, (H0_local - H0_CMB)/H0_CMB ~ 0.084."""
    return (H0_LOCAL - H0_CMB) / H0_CMB


def void_fraction_to_match_tension() -> float:
    """The void fraction at which the timescape offset R(f_v)-1 equals the
    observed tension fraction. Solve 4 f^2 + f + 4 = (1+tau) 2 (2+f) for f.
    Comes out ~0.47 -- BELOW the observed 0.76, so the real void fraction has
    headroom to spare."""
    tau = observed_tension_fraction()
    # 4 f^2 + (1 - 2(1+tau)) f + (4 - 4(1+tau)) = 0
    a = 4.0
    b = 1.0 - 2.0 * (1.0 + tau)
    c = 4.0 - 4.0 * (1.0 + tau)
    disc = b * b - 4.0 * a * c
    return (-b + disc**0.5) / (2.0 * a)


def resolves_in_right_direction(f_v0: float = FV0_BESTFIT) -> bool:
    """The mechanism resolves the tension iff the offset is positive (local >
    global) and at least as large as the observed gap."""
    return wall_vs_average_offset(f_v0) >= observed_tension_fraction() > 0.0


def early_universe_unaffected(tol: float = 1e-3) -> bool:
    """The clock effect is negligible at recombination (R-1 < tol there), so the
    CMB sound-horizon H0 is left intact -- the tension is a late-time artifact."""
    f_rec = void_fraction_at_recombination()
    return (dressed_over_bare(f_rec) - 1.0) < tol
