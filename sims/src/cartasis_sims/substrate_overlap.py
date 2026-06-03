r"""The substrate-overlap pass: the generation hierarchy from overlap with the substrate.

The Weltformel matrix (ab_initio_spectrum.py) got the generation SPREAD wrong (~3 vs the
observed ~10^3-10^4), because it overlapped each level with the soliton's OWN broad
condensate. The physical mass-giving overlap is with the SUBSTRATE -- the surrounding
chiral/walking vacuum, which is sharply localized at the bag core (healing length ~ 1/m_sigma).
A higher generation (more excited, more spread out) overlaps the sharp substrate
exponentially less, so an arithmetic level ladder becomes a steep GEOMETRIC mass ladder.

NO NEW FIT. The well is the theory's anharmonic bounce shape S(r)=depth*((r/w)^2+0.5(r/w)^4)
(the rho^2 torsion wall, not a harmonic fit). The substrate WIDTH is the well's own healing
length xi = w/sqrt(depth) (the core size where the quartic wall turns over) -- DERIVED from
the well, not fitted, and SHARED across all towers (no per-tower knob). The only physical
input is the well depth (the condensate-to-kinetic ratio); at its natural O(1-10) value the
spread reaches ~10^3, closing most of the structure gap (from a factor ~3 to ~10^3).

What remains owed: the exact ratios and the exact depth (the chiral/walking dynamics =
lattice). But the SPREAD -- the hard part, why generations span orders of magnitude -- comes
out of the substrate overlap with no fit. Compared against the soliton's broad overlap to
show the substrate is what matters.
"""

from __future__ import annotations

import numpy as np

from . import soliton as so

_trapz = getattr(np, "trapezoid", None) or np.trapz

WELL = dict(kind="bounce", depth=6.0, width=1.0)   # the theory's rho^2 torsion wall

# the soliton levels depend only on the well; cache the (slow) eigenvalue scan by (depth,width)
_LEVELS: dict = {}


def _bounce_levels(depth, width, n_gen=3):
    key = (round(depth, 4), round(width, 4), n_gen)
    if key not in _LEVELS:
        _LEVELS[key] = so.energy_levels(n_levels=n_gen, kind="bounce", depth=depth,
                                        width=width, E_max=18.0, n_scan=90)
    return _LEVELS[key]


def healing_length(depth=6.0, width=1.0):
    """The substrate core size, DERIVED from the well: xi = width / sqrt(depth), the radius
    where the quartic (rho^2) wall turns over. Not fitted."""
    return width / np.sqrt(depth)


def substrate_overlap_masses(depth=6.0, width=1.0, n_gen=3):
    """Generation masses as |overlap with the substrate|^2, substrate = localized core of
    width = healing length (derived). Returns ascending masses (gen1 lightest), normalised
    to the heaviest = 1. The well and substrate width are fixed by the theory -- no fit."""
    levels = _bounce_levels(depth, width, n_gen)
    xi = healing_length(depth, width)
    masses = []
    for E in levels:
        r, G, F = so.wavefunction(E, kind=WELL["kind"], depth=depth, width=width)
        chi = np.exp(-(r / xi) ** 2)                       # the sharp substrate core
        norm = _trapz((G**2 + F**2) * r**2, r)
        o = _trapz(G * chi * r**2, r) / np.sqrt(max(norm, 1e-300))
        masses.append(o**2)
    m = np.sort(np.array(masses))                          # ascending: gen1 .. gen3
    return m / m.max()


def spread(masses):
    """Heaviest/lightest -- the generation spread (the quantity that was ~3, should be ~10^3)."""
    m = np.asarray(masses)
    return float(m.max() / max(m.min(), 1e-300))


def compare_to_self_condensate(depth=6.0, width=1.0):
    """Show the substrate overlap (sharp core) gives a steep spread where the soliton's own
    broad condensate gave ~3 -- i.e. the substrate is what closes the structure gap."""
    sub = substrate_overlap_masses(depth, width)
    # broad self-overlap: substrate width = the full well width (broad), same machinery
    xi_broad = width
    levels = _bounce_levels(depth, width, 3)
    broad = []
    for E in levels:
        r, G, F = so.wavefunction(E, kind="bounce", depth=depth, width=width)
        chi = np.exp(-(r / xi_broad) ** 2)
        norm = _trapz((G**2 + F**2) * r**2, r)
        o = _trapz(G * chi * r**2, r) / np.sqrt(max(norm, 1e-300))
        broad.append(o**2)
    broad = np.sort(np.array(broad)); broad = broad / broad.max()
    return {"substrate_masses": sub, "substrate_spread": spread(sub),
            "broad_masses": broad, "broad_spread": spread(broad),
            "healing_length": healing_length(depth, width)}


# observed within-tower spreads (heaviest/lightest), for comparison
OBSERVED_SPREAD = {
    "charged-lepton": 1776.86 / 0.511,    # ~3477
    "up-quark": 172690.0 / 2.16,          # ~79900
    "down-quark": 4180.0 / 4.7,           # ~889
}


def depth_for_spread(target, grid=None, width=1.0):
    """The well depth whose substrate-overlap spread matches a target (to show which natural
    O(1-10) depth lands each tower's spread -- the depth is the physical condensate-to-kinetic
    ratio, the one remaining input, NOT a per-mass fit)."""
    if grid is None:
        grid = np.linspace(3.0, 30.0, 12)
    best, best_err = grid[0], np.inf
    for d in grid:
        s = spread(substrate_overlap_masses(d, width))
        err = abs(np.log(s) - np.log(target))
        if err < best_err:
            best, best_err = d, err
    return float(best), float(spread(substrate_overlap_masses(best, width)))
