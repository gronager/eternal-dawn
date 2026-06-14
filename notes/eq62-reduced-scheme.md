# Lab note — the reduced scheme for Eq 6.2: solve once, turn the A4 crank (2026-06-14)

Status: **lab notes.** The canonical statement of how the A4 symmetry collapses the three-generation
problem, tying the verified pieces together. Supersedes the scattered tetrahedral/Koide threads as the
*method*.

## The claim
The three charged leptons (and three neutrinos) are NOT three radial levels of Eq 6.2. They are the **one**
bound level of Eq 6.2 placed in the **three components of the A4 triplet 3**, degenerate by Schur until A4
breaks. So "the well binds one level" was never a failure — it is the singlet, exactly as it must be. The
threeness is the flavour index, not the radial quantum number. (This is why every radial tower and every
topological tower failed to give a clean 3: the threeness is not spatial.)

## Two symmetries, two reductions
- **Rotational (spatial):** 3D Dirac -> the 1D radial system (grand spin K). This is eq:radial.
- **A4 (flavour):** the THREE coupled generation equations -> ONE. By Schur (3 irreducible), all three obey
  the same radial equation with the same eigenvalue in the symmetric limit. Solve it once.

## The reduced program (the "5 pages" gone)
1. **One differential equation.** Solve the A4-singlet self-consistent soliton (Eq 6.2, first order, keeps
   the small component F):
       dG/dr = -(k/r)G + (1/hc)(E + Mc^2 + S)F
       dF/dr = +(k/r)F - (1/hc)(E - Mc^2 - S)G
   with the scalar self-energy S(r) = -g sigma(r) sourced self-consistently by the bound level (the chiral
   condensate the field digs), iterated to a fixed point. Output: the eigenvalue E0 and the wavefunction
   (G0, F0); the configurational-mass amplitude m0 = integral (G0^2 - F0^2) M(r) dr.
   [code: dirac_soliton.solve_dirac_soliton / self_consistent.solve_soliton / well_spectrum]
2. **Pure group theory (no calculus).** The three generations are psi_i = (G0,F0) x e_i. A4-breaking enters
   as a perturbation delta_M in the non-singlet channels (1' + 1'' + 3). First-order PT + the A4
   Clebsch-Gordan force the AMPLITUDE (sqrt-mass) matrix to be CIRCULANT, so
       sqrt(m_k) = m0 (1 + sqrt2 cos(delta + 2 pi k/3)),   m_k = (sqrt m_k)^2,   k=0,1,2.
   The sqrt2 and the 120 degrees are representation theory, NOT solved for.
   [code: a4_symmetry (proves circulant is forced), a4_flavor / koide (the pattern)]
3. **The magnitude.** The hierarchy is the NODE: 1 + sqrt2 cos(delta + ...) ~ 0 for the electron -- the
   near-critical amplification (criticality_scan), not a separate calculation.

## What is solved vs input
- **Solved:** the radial profile and m0 (one 1D self-consistent eigenvalue problem); the spectrum FORM
  (sqrt2 / 120 / circulant / TBM, all from A4).
- **Input (the residue):** the breaking phase delta (= the criticality phase, possibly CNS-selected) and
  the breaking scale; the small TBM-breaking theta13. Three masses from (m0, delta); form fixed by A4.

## The honest caveat
A4 is a flavour symmetry: it cuts the number of equations (3 -> 1), not the spatial dimension (that is the
rotational 1D reduction, separate). Between the one radial solve and the final spectrum everything is
algebra. The dynamical selection of A4 over U(3)/S4 is the four-body VEV-alignment question (four-fermion ->
tetrahedral cluster); the symmetry proof (a4_symmetry) is the other half.

## Module map
- Step 1 (radial solve): `dirac_soliton`, `self_consistent`, `well_spectrum`, `chiral_soliton`.
- Step 2 (A4 -> circulant/Koide + TBM): `a4_symmetry` (rigorous invariance), `a4_flavor`, `koide`.
- Step 3 (criticality magnitude): `condensate_excitation.criticality_scan`.
- Mixing readout: `mixing` (CKM aligned, PMNS tri-bimaximal, the 2+1 oscillation spectrum).
