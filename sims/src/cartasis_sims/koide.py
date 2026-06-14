r"""Koide's relation as the fingerprint of a Z3-symmetric three-sector (the topological generations).

Koide (1983): for the charged leptons  Q = (m_e+m_mu+m_tau)/(sqrt m_e+sqrt m_mu+sqrt m_tau)^2 = 2/3,
to 1e-5. Foot's geometric reading: the three sqrt-masses are three points 120 degrees apart on a
circle, amplitude sqrt2 about their mean,

    sqrt(m_k) = m0 * (1 + sqrt2 * cos(delta + 2 pi k/3)),   k = 0,1,2,

with the sqrt2 amplitude FIXED (that is exactly what makes Q = 2/3 -- see below), so just TWO numbers
(scale m0, phase delta) fix all three masses. Algebra: sum_k cos = 0 and sum_k cos^2 = 3/2, so
sum sqrt m = 3 m0 and sum m = m0^2(3 + 2*sqrt2*0 + 2*3/2) = 6 m0^2, giving Q = 6/9 = 2/3 for ANY delta.

Framework reading (the synthesis of Routes 1 & 3): three generations = three Z3-related topological
sectors (Hopfion charges sharing one condensate). The Z3 forces the sqrt(configurational-mass)
overlaps to sit at 120-degree phases -> Koide. The HIERARCHY (3477) is NOT a power law: it is the
phase delta placing one sector near the zero of (1 + sqrt2 cos), where the overlap nearly cancels ->
the anomalously tiny electron. Small delta-shifts move the (near-node) electron hugely while barely
touching mu, tau -- the steep sensitivity the configurational overlap supplies near a node.

Established: Koide (1983), Foot (1994, the 45-degree/120-degree geometry). New here: the identification
of the Z3 phases with distinct topological sectors, so Koide is the topological-generation fingerprint
and a 2-parameter mass law per tower.
"""
from __future__ import annotations

import numpy as np

# masses in MeV
LEPTONS = np.array([0.51099895, 105.6583755, 1776.86])
UP_QUARKS = np.array([2.16, 1270.0, 172_690.0])           # u, c, t (current, illustrative)
DOWN_QUARKS = np.array([4.67, 93.4, 4180.0])              # d, s, b


def koide_ratio(m):
    """Q = (sum m)/(sum sqrt m)^2. Charged leptons -> 2/3; bounded in [1/3, 1]."""
    m = np.asarray(m, float)
    s = np.sqrt(m)
    return float(m.sum() / s.sum() ** 2)


def z3_fit(m):
    """Fit sqrt(m_k)=m0(1+sqrt2 cos(delta+2pi k/3)) with the sqrt2 amplitude FIXED. m0=mean(sqrt m);
    scan delta and the k-assignment for the best phase. Returns (m0, delta, perm, predicted masses)."""
    from itertools import permutations
    m = np.asarray(m, float)
    sq = np.sqrt(m)
    m0 = sq.mean()
    best = None
    for perm in permutations(range(3)):
        for d in np.linspace(0, 2 * np.pi, 4000):
            pred = m0 * (1 + np.sqrt(2) * np.cos(d + 2 * np.pi * np.array(perm) / 3))
            err = np.sum((pred - sq) ** 2)
            if best is None or err < best[0]:
                best = (err, d, perm)
    _, d, perm = best
    pred = (m0 * (1 + np.sqrt(2) * np.cos(d + 2 * np.pi * np.array(perm) / 3))) ** 2
    return m0, float(d), perm, pred


def node_proximity(m0, delta, perm):
    """How close each sector sits to the overlap zero (1+sqrt2 cos)=0 -- the small value = the light
    generation. Returns the three (1+sqrt2 cos) factors; the near-zero one is the anomalously light rung."""
    k = np.array(perm)
    return 1 + np.sqrt(2) * np.cos(delta + 2 * np.pi * k / 3)


# neutrino mass-squared splittings (NuFIT/PDG 2024), eV^2
DM21 = 7.42e-5
DM31_NO = 2.515e-3       # normal ordering, m1 lightest
DM32_IO = -2.498e-3      # inverted ordering, m3 lightest


def neutrino_koide_range(ordering="NO", n=4000, ml_max=0.20):
    """Scan the lightest neutrino mass and return the achievable Koide-Q range given the measured
    splittings, the closest approach to 2/3, and the implied sum of masses. The point: with the
    measured Delta m^2, Q=2/3 is NOT reachable for neutrinos (NO tops at ~0.586, IO at ~0.500) -- the
    exact-Koide (sqrt2-amplitude Z3) does not extend to the neutral sector. This is CONSISTENT with
    neutrinos being the most Z3-broken sector (largest mixing, PMNS): exact (leptons, no mixing) ->
    approximate (quarks, CKM) -> unreachable (neutrinos). Returns dict with Q-range, best point, Sum m."""
    ml = np.linspace(0.0, ml_max, n)
    if ordering == "NO":
        m1, m2, m3 = ml, np.sqrt(ml**2 + DM21), np.sqrt(ml**2 + DM31_NO)
    else:
        m3 = ml
        m1 = np.sqrt(ml**2 + abs(DM32_IO) - DM21)
        m2 = np.sqrt(ml**2 + abs(DM32_IO))
    M = np.stack([m1, m2, m3], axis=0)
    Q = M.sum(0) / np.sqrt(M).sum(0) ** 2
    i = int(np.argmax(Q))
    return {"ordering": ordering, "Q_min": float(Q.min()), "Q_max": float(Q.max()),
            "reaches_2_3": bool(Q.max() >= 2 / 3), "ml_at_maxQ": float(ml[i]),
            "sum_m_at_maxQ": float(M[:, i].sum()), "masses_at_maxQ": M[:, i].copy()}


def report():
    print("Koide's relation as the Z3 three-sector fingerprint\n")
    for name, m in [("charged leptons", LEPTONS), ("up quarks", UP_QUARKS),
                    ("down quarks", DOWN_QUARKS)]:
        Q = koide_ratio(m)
        m0, d, perm, pred = z3_fit(m)
        fac = node_proximity(m0, d, perm)
        print(f"  {name:>16}:  Koide Q = {Q:.5f}  (2/3={2/3:.5f})   "
              f"delta={np.degrees(d):.1f} deg")
        print(f"    {'':16}   2-param fit (MeV): {np.round(np.sort(pred),3)}")
        print(f"    {'':16}   observed   (MeV): {np.round(np.sort(m),3)}")
        print(f"    {'':16}   (1+sqrt2 cos) factors = {np.round(np.sort(fac),3)}  "
              f"-> min={np.min(np.abs(fac)):.3f} = the light rung (overlap near-node)\n")
    print("  Charged leptons: Q=2/3 to 1e-5; the 2-param (scale, phase) law reproduces all three to")
    print("  <1%. The hierarchy lives in the PHASE delta (electron near the overlap node), not a power")
    print("  law. Framework: the Z3 = three topological sectors; the near-node = the configurational")
    print("  overlap nearly cancelling for the lightest generation.\n")

    print("  Neutrinos -- can exact Koide 2/3 hold given the measured Delta m^2?")
    for ordering in ("NO", "IO"):
        r = neutrino_koide_range(ordering)
        verdict = "REACHABLE" if r["reaches_2_3"] else "NOT reachable"
        print(f"    {ordering}: Q in [{r['Q_min']:.3f}, {r['Q_max']:.3f}] -> 2/3 {verdict}; "
              f"closest at m_light={r['ml_at_maxQ']*1e3:.1f} meV, Sum m={r['sum_m_at_maxQ']*1e3:.0f} meV")
    print("    => exact Koide does NOT extend to neutrinos -- CONSISTENT with the neutral sector being")
    print("    the most Z3-broken (largest mixing). Closest-to-Z3 = NORMAL ordering, hierarchical,")
    print("    Sum m ~ 59 meV (falsifiable: cosmology Sum m<120 meV, CMB-S4 ~20 meV sensitivity).")


if __name__ == "__main__":
    report()
