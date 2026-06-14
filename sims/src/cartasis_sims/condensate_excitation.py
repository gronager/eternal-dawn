r"""Generations as excitations of the CONDENSATE, not the bound state -- the ceiling-free hierarchy.

Every "bound-state excitation" picture (radial tower, winding zero-modes, Hopfion energy) capped the
mass span at ~2-15: three wavefunctions in ONE well overlap similarly. This module tests the other
idea -- the three generations are three EXCITATIONS of the condensate sigma(r), and the (fixed) fermion
sits in each. The mass is still the configurational overlap

    m_n = integral rho_f(r) sigma_n(r) d^3r ,   sigma_n = the n-th condensate radial mode (n nodes),

but now the well DIFFERS per generation. The decisive feature: as the fermion density rho_f matches
the condensate ground mode, its overlap with the EXCITED modes -> 0 by ORTHOGONALITY. So the span is
NOT capped -- it diverges at the orthogonality point -- and the near-cancellation IS Koide's node (the
anomalously light generation). One control parameter a/b (fermion width / condensate-mode width), like
the Koide phase.

Findings (report):
  * the span is ceiling-free: it diverges near a/b=1 (orthogonality) -- the first mechanism that can
    reach the lepton span 3477 naturally (bound-state towers cannot);
  * the lightest generation = the most-excited condensate mode (most cancellation) -> a clean rung
    assignment (electron = highest condensate excitation, tau = condensate ground);
  * the overlap masses PASS THROUGH Koide Q=2/3 at special a/b (compatible, not forced);
  * one parameter fits one ratio (the electron) but the muon lands ~3x light -- the specific oscillator
    modes are not the lepton pattern by themselves; the Z3/Koide structure is still needed to fix BOTH
    ratios. So this supplies the MAGNITUDE MECHANISM (steep, ceiling-free) that the bound-state pictures
    lacked; the Z3 supplies the pattern.
"""
from __future__ import annotations

import numpy as np
from scipy.special import genlaguerre

_trapz = getattr(np, "trapezoid", None) or np.trapz


def _grid(rmax=14.0, N=6000):
    return np.linspace(1e-4, rmax, N)


def cond_mode(n, b, r):
    """The n-th condensate radial mode (3D isotropic, s-wave; n radial nodes), width b, L2-normalised
    with the r^2 measure: sigma_n(r) ~ L_n^{1/2}(r^2/b^2) exp(-r^2/2b^2)."""
    x = r**2 / b**2
    R = genlaguerre(n, 0.5)(x) * np.exp(-x / 2)
    return R / np.sqrt(_trapz(R**2 * r**2, r))


def overlap_masses(a, b=1.0, nmax=3, r=None):
    """Configurational overlaps m_n = integral rho_f sigma_n r^2 dr for a fixed fermion density
    rho_f ~ exp(-r^2/2a^2) and condensate modes n=0..nmax-1. Returns the signed overlaps."""
    if r is None:
        r = _grid()
    rho = np.exp(-r**2 / (2 * a**2))
    rho /= _trapz(rho * r**2, r)
    return np.array([_trapz(rho * cond_mode(n, b, r) * r**2, r) for n in range(nmax)])


def span(a, b=1.0, r=None):
    """Mass span max/min over the three generations (|overlap|)."""
    m = np.abs(overlap_masses(a, b, 3, r))
    return float(m.max() / max(m.min(), 1e-300))


def koide_Q(a, b=1.0, square=False, r=None):
    """Koide ratio of the condensate-overlap masses (mass = |overlap|, or overlap^2 if square)."""
    m = np.abs(overlap_masses(a, b, 3, r))
    if square:
        m = m**2
    s = np.sqrt(m)
    return float(m.sum() / s.sum() ** 2)


def report():
    print("Generations as condensate excitations -- the ceiling-free hierarchy\n")
    r = _grid()
    print("  a/b    span(0..2)     Koide-Q(|m|)   m0:m1:m2 (|overlap|, normalised)")
    for ab in (0.6, 0.75, 0.85, 0.9, 0.95, 1.05, 1.15, 1.3):
        m = np.abs(overlap_masses(ab, 1.0, 3, r))
        mn = m / m.max()
        print(f"  {ab:>4.2f}   {span(ab,1.0,r):>10.1f}    {koide_Q(ab,r=r):>8.4f}      "
              f"{mn[0]:.4f} : {mn[1]:.4f} : {mn[2]:.5f}")
    print("\n  - span DIVERGES near a/b=1 (fermion orthogonal to excited modes) -> no ceiling, reaches")
    print("    3477 naturally; the near-cancellation = Koide's node = the anomalously light generation.")
    print("  - Koide Q crosses 2/3 at special a/b (compatible, not forced).")
    print("  - lightest gen = most-excited condensate mode (most cancellation): electron = highest")
    print("    condensate excitation, tau = condensate ground. One knob a/b fits one ratio; the Z3/Koide")
    print("    structure fixes both. This supplies the MAGNITUDE mechanism the bound-state towers lacked.")


if __name__ == "__main__":
    report()
