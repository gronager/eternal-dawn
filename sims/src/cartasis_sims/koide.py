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
    print("  overlap nearly cancelling for the lightest generation.")


if __name__ == "__main__":
    report()
