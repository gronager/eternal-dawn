r"""What the transition diagrams (CKM/PMNS) and neutrino oscillations tell us about the knot field.

The charged-current vertices are knot-reconnection events. Reading all three sectors together:

 - CHARGED LEPTONS: the W only does l <-> nu_l (within a generation); there is NO cross-generation
   charged-lepton vertex (no mu->e gamma). tau decay is tau -> nu_tau + W, W -> mu nu_mubar: the tau
   becomes its OWN nu_tau and the W independently makes the lighter pair -- no direct tau->mu vertex.
   The charged-lepton topology is RIGID.
 - QUARKS: cross-generation vertices EXIST but are CKM-suppressed (s->u, b->c, b->u ...). The W can
   re-tie a quark knot across generations. The quark topology is SOFT (weakly reconnectable).
 - NEUTRINOS: the mass eigenstates (the true knots) are misaligned with the flavor (W-coupling) basis
   by the LARGE PMNS angles, so they oscillate.

The unifying statement -- MIXING = MISALIGNMENT of the mass basis (true knots) and the W-coupling
(flavor) basis, set by ANCHORING:
  * charged knots (up, down, charged leptons) are anchored by charge -> mass & flavor bases ALIGN
    -> small mixing (CKM ~ identity);
  * the neutral knot (neutrino) is unanchored -> its orientation FLOATS -> strong misalignment
    -> large mixing (PMNS ~ tri-bimaximal). The same charge-anchoring that gives the charged leptons
    EXACT Koide 2/3 (no mixing) lets the neutral knots float (Koide fails, large mixing).

Neutrino oscillation = a literal INTERFEROMETER of the knot energies (phase exp(-i m_i^2 L/2E)):
  * TWO frequencies (Dm21^2, Dm31^2) = the two gaps of a 3-knot spectrum -> 3 knots, no light sterile;
  * a 2+1 spectrum (Dm31^2/Dm21^2 ~ 34): two near-degenerate knots + one split off;
  * near-maximal theta_23 ~ a mu-tau exchange symmetry of the knot field;
  * near-tri-bimaximal PMNS = the A4/Z3 discrete symmetry of the three knots -- the SAME Z3 that Koide
    reads off the masses. Oscillations measure the OFF-diagonal (mixing); Koide the diagonal (masses);
    together they pin the three-knot Z3/A4 structure, and both favour NORMAL ordering + hierarchical
    Sum m ~ 60 meV.
"""
from __future__ import annotations

import numpy as np

# oscillation gaps (NuFIT/PDG 2024), eV^2
DM21 = 7.42e-5
DM31 = 2.515e-3
# best-fit mixing angles (deg): PMNS (leptons) and CKM (quarks)
PMNS_DEG = (33.4, 49.0, 8.6)     # theta12, theta23, theta13
CKM_DEG = (13.0, 2.4, 0.20)


def mass_gap_structure():
    """The 'double oscillation' as the knot spectrum: 2 gaps = 3 knots; the 2+1 (two close, one far)."""
    return {"Dm21": DM21, "Dm31": DM31, "ratio": DM31 / DM21,
            "n_knots": 3, "n_frequencies": 2, "pattern": "2+1 (two near-degenerate + one split)"}


def mixing_matrix(t12_deg, t23_deg, t13_deg):
    """Standard CKM/PMNS parametrisation |U| from three angles (CP phase set to 0)."""
    t12, t23, t13 = np.radians([t12_deg, t23_deg, t13_deg])
    s, c = np.sin, np.cos
    R23 = np.array([[1, 0, 0], [0, c(t23), s(t23)], [0, -s(t23), c(t23)]])
    R13 = np.array([[c(t13), 0, s(t13)], [0, 1, 0], [-s(t13), 0, c(t13)]])
    R12 = np.array([[c(t12), s(t12), 0], [-s(t12), c(t12), 0], [0, 0, 1]])
    return R23 @ R13 @ R12


def tribimaximal():
    """The TBM matrix: sin^2 theta12 = 1/3, theta23 = 45, theta13 = 0 (the A4/Z3 ideal)."""
    return mixing_matrix(np.degrees(np.arcsin(1 / np.sqrt(3))), 45.0, 0.0)


def offdiagonal_power(U):
    """A crude misalignment measure: 1 - mean|diagonal| of |U| (0 = aligned/identity, large = mixed)."""
    return float(1 - np.mean(np.abs(np.diag(np.abs(U)))))


def report():
    g = mass_gap_structure()
    print("Transitions & oscillations -- the tell about the knot field\n")
    print(f"  the 'double oscillation' = {g['n_frequencies']} gaps of {g['n_knots']} knots: "
          f"Dm31/Dm21 = {g['ratio']:.0f} -> {g['pattern']}")
    print("  (a 4th light knot would add a 3rd frequency -> none seen -> no light sterile)\n")
    pmns = mixing_matrix(*PMNS_DEG); ckm = mixing_matrix(*CKM_DEG); tbm = tribimaximal()
    print(f"  PMNS off-diagonal power = {offdiagonal_power(pmns):.2f}  (large; |PMNS-TBM|max="
          f"{np.abs(np.abs(pmns)-np.abs(tbm)).max():.2f} -> near tri-bimaximal = A4/Z3)")
    print(f"  CKM  off-diagonal power = {offdiagonal_power(ckm):.3f}  (small; ~ identity = aligned)")
    print("\n  MIXING = MISALIGNMENT of mass (true-knot) vs flavour (W-coupling) basis, set by ANCHORING:")
    print("    charged knots anchored by charge -> aligned -> small mixing (CKM, and Koide exact);")
    print("    neutral knot (neutrino) unanchored -> floats -> large mixing (PMNS), Koide fails.")
    print("  Oscillations read the OFF-diagonal (Z3 via mixing); Koide the diagonal (Z3 via masses);")
    print("  both pin the three-knot Z3/A4 and favour normal ordering + hierarchical Sum m ~ 60 meV.")


if __name__ == "__main__":
    report()
