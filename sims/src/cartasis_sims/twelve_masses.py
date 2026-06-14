r"""All 12 fermion masses as 4 towers x the A4 circulant crank -- and an honest parameter count.

Each tower (charged leptons, up quarks, down quarks, neutrinos) is a generation triplet whose
sqrt-mass matrix is the A4 Z3-circulant, eigenvalues

    sqrt(m_k) = c0 (1 + A cos(delta + 2 pi k/3)),   k = 0,1,2,

with c0 the scale, delta the phase, A the amplitude. A4 (unbroken) FORCES A = sqrt2 (= Koide Q=2/3).

Findings (report):
  * charged leptons: A = sqrt2 EXACTLY (Koide 2/3) -> A4-exact -> 3 masses from 2 numbers (c0, delta);
  * up/down quarks and neutrinos: A deviates from sqrt2 (+24% / +9% / -14%) -> A4 BROKEN (the same
    breaking that gives CKM/PMNS mixing) -> 3 params each, no net mass prediction (but the circulant
    STRUCTURE still organises them, and TBM mixing is predicted for the neutral sector).

Parameter accounting (masses only):
  * IF A4 were exact in every tower:  4 x (c0, delta) = 8 params for 12 masses  (12 -> 8, the goal);
  * AS IT STANDS (only leptons exact): 2 + 3 + 3 + 3 = 11 params (12 -> 11, the lepton Koide the one
    clean mass prediction). The 3 'lost' predictions are the quark/neutrino A4-breakings -- tied to the
    mixing, so relating A-deviation to CKM/PMNS would recover toward 12 -> 8. The BIG reduction is the
    MIXING: A4 predicts tri-bimaximal (PMNS) ~ for free vs the SM's ~8 mixing inputs.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

SQRT2 = np.sqrt(2.0)

# masses in MeV; neutrinos: normal-ordering, hierarchical (m1~0), in meV scaled to MeV here
TOWERS = {
    "charged leptons": np.array([0.51099895, 105.6583755, 1776.86]),
    "up quarks": np.array([2.16, 1270.0, 172690.0]),
    "down quarks": np.array([4.67, 93.4, 4180.0]),
    "neutrinos (NO)": np.array([1e-6, 8.6e-3, 5.0e-2]),   # ~meV, m1->0 (only Delta m^2 known)
}


def fit_circulant(masses):
    """Fit (c0, amplitude A, delta) of sqrt(m_k)=c0(1+A cos(delta+2pi k/3)) to the three masses.
    Returns (c0, A, delta_deg, Koide_Q). A4-exact <=> A = sqrt2."""
    sq = np.sort(np.sqrt(np.asarray(masses, float)))
    c0 = sq.mean()

    def resid(p):
        c0, A, d = p
        pred = np.sort([c0 * (1 + A * np.cos(d + 2 * np.pi * k / 3)) for k in range(3)])
        return pred - sq

    best = None
    for d0 in np.linspace(0, 2 * np.pi, 36):
        r = least_squares(resid, [c0, 1.0, d0])
        if best is None or r.cost < best.cost:
            best = r
    c0, A, d = best.x
    Q = float(np.sum(masses) / np.sum(np.sqrt(masses)) ** 2)
    return float(c0), float(abs(A)), float(np.degrees(d) % 360), Q


def report():
    print("All 12 fermion masses = 4 towers x the A4 circulant\n")
    print(f"  {'tower':16} {'scale c0':>10} {'amplitude A':>12} {'A/sqrt2':>9} {'Koide Q':>9}  status")
    n_params = 0
    for name, m in TOWERS.items():
        c0, A, d, Q = fit_circulant(m)
        ratio = A / SQRT2
        exact = abs(A - SQRT2) < 0.05
        n_params += 2 if exact else 3
        status = "A4-EXACT (predicts)" if exact else f"broken {(ratio-1)*100:+.0f}% (fit)"
        print(f"  {name:16} {c0:>10.3f} {A:>12.3f} {ratio:>9.3f} {Q:>9.4f}  {status}")
    print(f"\n  A4 forces A=sqrt2 (Koide 2/3) for an UNBROKEN tower -> 3 masses from (c0, delta).")
    print(f"  Mass params now: {n_params} for 12 masses (12 -> {n_params}); only the LEPTON tower is")
    print(f"  A4-exact (the one clean mass prediction). If A4 were exact in all 4: 8 params (12 -> 8).")
    print(f"  The quark/neutrino A-deviations ARE the A4-breaking = the CKM/PMNS mixing; relating them")
    print(f"  would recover toward 12 -> 8. The big win is the MIXING: A4 predicts tri-bimaximal PMNS.")


if __name__ == "__main__":
    report()
