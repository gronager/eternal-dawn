r"""Generations: an arithmetic ladder of internal levels -> a geometric mass hierarchy.

The three generations (e/mu/tau, and the quark towers) are the program's quarantined
gap: no ordinary ladder gives the observed hierarchy, which spans orders of magnitude.
This module shows the mechanism that *could* bridge it -- "Yukawa as overlap."

THE IDEA. Suppose the generations are internal excitation levels of the same soliton
(same colour, charge, weak isospin; different internal rung). The level ENERGIES are
roughly arithmetic (evenly spaced, E_n ~ E_0 + n*dE). But the observable MASS is not
the level energy -- in the soliton picture it is a wavefunction OVERLAP with the
localized mass-giving region (the condensate / the would-be Higgs), and overlaps are
EXPONENTIALLY sensitive to the level. A higher rung sits further out / penetrates a
barrier less, so its overlap with the localized core falls like

    O_n  ~  exp(-c n),      mass_n ~ |O_n|^2 ~ exp(-2 c n).

So an ARITHMETIC ladder in n becomes a GEOMETRIC ladder in mass: ln(mass) linear in n,
each generation a fixed multiple of the next. Equivalently: the "unnaturally tiny"
electron Yukawa (~3e-6) is just exp(-2c*0_offset) -- the exponential of an order-one
overlap suppression, not a fine-tuned small number. That is the qualitative win the
Standard Model lacks: a large hierarchy from O(1) inputs.

WHAT IS AND ISN'T DELIVERED. Computed (here, from the real soliton wavefunctions):
the overlap falls exponentially, producing a geometric hierarchy whose steepness is
set by the source size. Observed: the real fermion masses ARE approximately evenly
spaced in log (geometric-ish), exactly the structure this produces. NOT delivered:
the precise ratios (leptons 207, 17 -- a decreasing, not constant, ratio), which need
the detailed level/overlap structure and are the owed refinement. The mechanism
explains the SHAPE of the hierarchy, not yet its exact numbers.
"""

from __future__ import annotations

import numpy as np

from . import soliton as so

_trapz = getattr(np, "trapezoid", None) or np.trapz

# Standard-Model charged-fermion masses (MeV), PDG-ish, for the log-spacing comparison.
SM_FERMIONS = {
    "leptons (e, mu, tau)": [0.511, 105.66, 1776.86],
    "up-type (u, c, t)": [2.16, 1270.0, 172690.0],
    "down-type (d, s, b)": [4.67, 93.4, 4180.0],
}


def core_overlap(E, source_size=0.25, kind="linear", depth=6.0, width=1.0):
    """Overlap of the level-E internal wavefunction with a localized core source
    chi(r) = exp(-(r/source_size)^2): O = |<chi|G>| / ||psi||."""
    r, G, F = so.wavefunction(E, kind=kind, depth=depth, width=width)
    chi = np.exp(-(r / source_size) ** 2)
    norm = _trapz((G**2 + F**2) * r**2, r)
    o = _trapz(G * chi * r**2, r)
    return abs(o) / np.sqrt(max(norm, 1e-300))


def overlap_masses(levels, source_size=0.25, kind="linear", depth=6.0, width=1.0):
    """mass_n ~ |O_n|^2 for a set of internal levels (the overlap hierarchy)."""
    ov = np.array([core_overlap(E, source_size, kind, depth, width) for E in levels])
    return ov**2


def self_consistent_hierarchy(g=14.0, m_sigma=0.5, ngen=3, substrate=0.4, R=14.0, N=1200,
                              iters=900, mix=0.12):
    r"""The generation hierarchy from the REAL self-consistent torsiton well (not the toy).

    Solve the self-consistent four-fermion bag (self_consistent.solve_soliton), read off the lowest
    `ngen` radial rungs (the generations: same charge/colour/isospin, different radial excitation),
    and form the configurational mass of each as a DENSITY overlap m_n = \int u_n^2 chi(r) dr -- the
    structure of Eq.~(configmass). Two substrates are returned: the BROAD condensate sigma(r), and a
    SHARP core of width `substrate` (the chiral-restored core, the localized mass-giving region of
    Section 11.5). The decisive result: the level energies E_n are ~arithmetic, the broad-overlap
    masses span only ~3, but the sharp-overlap masses span ~10^2-10^3 -- the arithmetic level ladder
    becomes the geometric MASS ladder, the observed size of the generation hierarchy, with NO fit to
    any fermion mass. The sharp spread SATURATES as the substrate narrows, so it is a robust limit.

    Returns a dict: levels E_n, the broad and sharp overlap masses (normalised, lightest=1), their
    spreads, and the number of bound rungs found."""
    from . import self_consistent as sc
    from .self_consistent import _solve_levels

    out = sc.solve_soliton(m0=1.0, g=g, m_sigma=m_sigma, n_fermions=1, n_levels=max(10, ngen + 2),
                           R=R, N=N, iters=iters, mix=mix)
    r = out["r"]
    h = r[1] - r[0]
    E, u = _solve_levels(out["M"], r, h, max(10, ngen + 2))
    n_bound = int(np.sum(E < 1.0))                    # rungs below the constituent mass m0=1
    n = min(ngen, max(n_bound, 1))

    def density_overlap(chi):
        m = np.array([_trapz(u[:, k] ** 2 * chi, r) for k in range(n)])
        m = np.abs(m)
        return m / m.min()

    broad = density_overlap(out["sigma"])
    sharp = density_overlap(np.exp(-(r / substrate) ** 2))
    return {
        "levels": np.asarray(E[:n]),
        "level_spacings": np.diff(E[:n]),
        "broad_masses": broad,
        "sharp_masses": sharp,
        "broad_spread": float(broad.max() / broad.min()),
        "sharp_spread": float(sharp.max() / sharp.min()),
        "n_bound": n_bound,
        "converged": bool(out["converged"]),
    }


def geometric_factor(masses):
    """The per-generation mass ratio from a log-linear fit: each rung is `factor` x
    the next. A constant factor == a clean geometric hierarchy."""
    n = np.arange(len(masses))
    slope = np.polyfit(n, np.log(masses), 1)[0]
    return float(np.exp(-slope))


def log_spacings(masses):
    """Successive ln(mass) gaps. Equal gaps == geometric. The SM fermions have gaps
    that are equal to within a factor ~2 -- approximately geometric."""
    lm = np.log(np.asarray(masses, dtype=float))
    return np.diff(lm)


def is_approximately_geometric(masses, tol=2.5):
    """True if the ln-mass gaps are equal to within a factor `tol` (geometric-ish)."""
    g = np.abs(log_spacings(masses))
    return float(np.max(g) / np.min(g)) < tol
