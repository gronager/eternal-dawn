# Gate calculation — the gravity–torsion soliton (Parts II/III)

First-approximation results for the load-bearing question of Parts II–III: does the
gravity–torsion potential bind, what is its spectrum, and can the electroweak
S parameter survive? Code: `sims/src/cartasis_sims/soliton.py`; figure:
`figures/scripts/draft_soliton_spectrum.py` → `figures/pdf/soliton_spectrum.pdf`.
SymPy-verified: the radial Dirac reduction, and the S-parameter sum-rule algebra.

## The framing the computation forced

The "potential" is **not external**. Integrating torsion out of Einstein–Cartan
leaves a Dirac field with a four-fermion (Hehl–Datta) self-interaction
`(ψ̄ γ5 γμ ψ)²`. So the well is the one the fermion's own axial-spin density digs —
the honest problem is a **nonlinear, self-consistent Dirac equation** (NJL / Soler
class), not a fixed 1-D well. The calculation here solves the radial Dirac equation
in a *fixed* well (first iteration); self-consistent back-reaction + Pauli filling is
the owed refinement and is where it gets harder.

## What came out

1. **Derrick is clearable.** A scalar lump is forbidden in 3-D (Derrick), but a Dirac
   (Soler) field evades it, and the two-terms-balance (attraction vs ρ²-torsion wall)
   is the same one that makes the cosmological bounce. So the soliton can exist; the
   explicit demonstration for the Part I potential is target #1.

2. **The spectrum is ~Regge.** A relativistic fermion in a confining scalar well gives
   a rising tower with `E²` linear in level index — to <1% for a linear (flux-tube)
   wall, bending up for the steeper anharmonic bounce wall. Matches the hadron Regge
   form, parameter-free (shape only).

3. **Scale caveat.** The slope is a *shape* in dimensionless units. The bare EC
   coupling (~G/c⁴) would put the absolute scale near the Planck mass; whether the
   particle-scale binding density equals the cosmological ρ_C is the open density
   question (target #5).

4. **S parameter — the make-or-break, NOT delivered.** Weinberg sum rules give
   `S = 4π (f_π/M_V)² (1 + M_V²/M_A²)`. QCD-like → **S ≈ 0.25** (the graveyard,
   above the LEP bound ≈ 0.1). S < 0.1 needs `M_V/f_π ≳ 13` (vs QCD's ≈ 8): a factor
   ∼1.6 of "walking"/near-conformal enhancement. The anharmonic bounce well *might*
   walk, but we have not computed the composite vector/axial spectrum, so S < 0.1 is
   **owed**, exactly as Part III's own status says.

## The tunnelling spanner (closes the loop to Part I)

The bounce/soliton well is finite → tunnelling both ways. **Outward** = decay /
Hawking evaporation (death). **Inward** = penetration to densities *above* ρ_C, which
is a region that bounces → a new interior, i.e. **baby-universe nucleation** (birth) —
the same instanton Part I uses (Ch. 3), now read as WKB through the membrane. The
membrane is a two-way barrier; the rates are wildly asymmetric (outward single-quantum,
inward instanton-suppressed ~10^{84}), which is why matter looks stable. In the
homogeneous cosmology the inward wall rises without bound (only quantum fuzz); it is
*local* collapse (soliton / black hole) that turns inward tunnelling into a real new
region.

## Honest verdict

The gate is **half-open**: the soliton plausibly exists (Derrick clearable) and its
tower is Regge (good). But the two numbers that would make the program compelling —
the absolute mass scale and S < 0.1 — are both still owed, and S in particular is a
real coin-flip the framework must compute and survive. Nothing here promotes Parts
II/III above "well-posed research program"; it sharpens the targets and confirms the
walls are where the drafts said.

## Update: degeneracy pressure (Pauli) — which way it pushes

`fermi_ball.py` (Thomas-Fermi, SymPy-verified Fermi-gas energy). Adding the Pauli
degeneracy KE to the self-bound drop (e(n) = e_kin - a n + b n^2):

- saturation density drops to ~81% and binding per particle to ~25% of the no-Pauli
  case — the drop gets **larger, less dense, and flatter-topped**.
- A flatter, more uniform (scale-flat) interior is the **walking** regime that
  *suppresses* the electroweak S. So degeneracy moves S in the FAVOURABLE direction.

It does NOT deliver S < 0.1 by itself (that needs the full self-consistent, colored,
composite V/A spectrum), but it confirms the intuition: Pauli is in the column that
helps. Figure: `figures/pdf/degeneracy.pdf`.
