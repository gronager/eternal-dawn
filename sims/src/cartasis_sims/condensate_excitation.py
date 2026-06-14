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


def real_soliton_modes(g=5.0, lam=6.0, nf=3, v=1.0, R=14.0, N=700, eps_break=0.4):
    """The HONEST test: the REAL chiral-soliton condensate modes, not toy oscillators. Solve the
    chiral soliton (chiral_soliton.solve_chiral), build the condensate-fluctuation operator
    H_sig = -d^2/dr^2 + lam(3 sigma_0^2 - v^2) (binding -lam v^2 in the bag, m_sigma^2=2 lam v^2 outside),
    find its BOUND modes (omega^2 < m_sigma^2), and overlap each with the fermion density. Returns
    (n_bound, overlaps, span, koide_Q). Finding: the real soliton binds ~2 modes (not 3); tuned to 3,
    they are near-DEGENERATE in overlap (span ~1.2, Q~1/3) -- the steep toy span needed the fermion to
    sit at near-orthogonality (a/b->1), a near-critical point the self-consistent soliton does NOT
    naturally occupy. So the mechanism is real but the magnitude wants a criticality driver, not free."""
    from scipy.linalg import eigh_tridiagonal
    from . import chiral_soliton as cs
    o = cs.solve_chiral(v=v, g=g, lam=lam, n_fermions=nf, R=R, N=N, eps_break=eps_break)
    r = o["r"]; h = r[1] - r[0]; s0 = o["sigma"]; rho = o["density"]
    thr = 2 * lam * v**2
    diag = 2.0 / h**2 + lam * (3 * s0**2 - v**2)
    off = -np.ones(len(r) - 1) / h**2
    w, U = eigh_tridiagonal(diag, off)
    ov = []
    for k in range(len(w)):
        if w[k] < thr * 0.999:
            u = U[:, k] / np.sqrt(h)
            ov.append(_trapz(rho * (u / r) * r**2, r))
    m = np.abs(np.array(ov))
    span_ = float(m.max() / m.min()) if len(m) else float("nan")
    Q = float(m.sum() / np.sqrt(m).sum() ** 2) if len(m) else float("nan")
    return len(ov), np.array(ov), span_, Q


def criticality_scan(g_values=None, lam=8.0, nf=2, eps_break=0.4, v=1.0, R=16.0, N=800, iters=400):
    """The criticality-driver test. Scan the coupling g toward marginal binding (the chiral phase
    boundary) and track the condensate-overlap span. Hypothesis: approaching the critical point drives
    the fermion to near-orthogonality with the excited condensate modes (the overlap NODE), so the span
    diverges. Finding: the span DOES grow sharply toward criticality (the driver is real and in the
    right direction), but the simple chiral soliton dissolves via a FIRST-ORDER transition (core flips
    abruptly) and caps the span at ~30 before the divergence. A SECOND-ORDER (continuous) chiral
    transition would let the system sit arbitrarily close to a/b->1 and deliver the full hierarchy.
    So the magnitude reduces to: IS THE CHIRAL TRANSITION (in this sector) FIRST OR SECOND ORDER?
    -- a real, lattice-testable question. Returns a list of dicts (g, core, margin, n_modes, span)."""
    from scipy.linalg import eigh_tridiagonal
    from . import chiral_soliton as cs
    if g_values is None:
        g_values = np.linspace(5.0, 2.64, 14)
    out = []
    for g in g_values:
        o = cs.solve_chiral(v=v, g=float(g), lam=lam, n_fermions=nf, R=R, N=N,
                            eps_break=eps_break, iters=iters)
        r = o["r"]; h = r[1] - r[0]; s0 = o["sigma"]; rho = o["density"]
        Mvac = abs(o["M"][-1]); E = np.array(o["E"]); Eb = E[(E > 0) & (E < Mvac)]
        margin = float(Eb.min() / Mvac) if len(Eb) else float("nan")
        diag = 2.0 / h**2 + lam * (3 * s0**2 - v**2); off = -np.ones(len(r) - 1) / h**2
        w, U = eigh_tridiagonal(diag, off)
        ov = [_trapz(rho * (U[:, k] / np.sqrt(h) / r) * r**2, r)
              for k in range(len(w)) if w[k] < 2 * lam * v**2 * 0.999]
        m = np.abs(np.array(ov))
        out.append({"g": float(g), "core": float(o["core_sigma"]), "margin": margin,
                    "n_modes": len(m), "span": float(m.max() / m.min()) if len(m) >= 2 else float("nan")})
    return out


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
    print("    structure fixes both. This supplies the MAGNITUDE mechanism the bound-state towers lacked.\n")

    print("  REAL chiral-soliton condensate modes (not toy oscillators):")
    n, ov, sp, Q = real_soliton_modes(g=5.0, lam=6.0, nf=3)
    print(f"    bound modes = {n}, overlaps = {np.round(ov,4)}, span = {sp:.1f}, Koide Q = {Q:.3f}")
    print("    -> the real soliton binds ~2-3 modes, and when 3 they are near-DEGENERATE (span ~1.2,")
    print("    Q~1/3), NOT steep. The toy's ceiling-free span needed a/b->1 (near-orthogonality) -- a")
    print("    near-critical point the self-consistent soliton does not naturally sit at. Mechanism real;")
    print("    magnitude wants a criticality driver (what pins a/b->1?), else it is the one tuned number.\n")

    print("  CRITICALITY-DRIVER scan (coupling g -> marginal binding / chiral phase boundary):")
    sc = criticality_scan(g_values=np.array([5.0, 4.0, 3.0, 2.8, 2.74, 2.68, 2.66]))
    print("    g     core_sig  margin(E/Mvac)  span")
    for d in sc:
        print(f"    {d['g']:.2f}   {d['core']:+.3f}     {d['margin']:.3f}        {d['span']:7.1f}")
    print("    -> span GROWS toward criticality (the driver is real); the excited-mode overlap heads to")
    print("    the NODE. But the chiral soliton dissolves FIRST-ORDER (~g=2.65), capping span ~30. A")
    print("    SECOND-ORDER transition would diverge -> magnitude = the ORDER of the chiral transition.")


if __name__ == "__main__":
    report()
