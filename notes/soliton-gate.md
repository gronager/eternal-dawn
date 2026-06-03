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

## Generations: arithmetic ladder -> geometric hierarchy (overlap mechanism)

`generations.py` + `figures/pdf/generations.pdf`. If generations are internal
excitation levels of the same soliton (same charge/colour/isospin, different rung),
the level energies are ~arithmetic but the observable mass is a wavefunction OVERLAP
with the localized condensate -- and overlaps are exponentially sensitive to the
level. Computed from the real soliton wavefunctions: the core overlap O_n falls
exponentially with n, so mass ~ |O_n|^2 is geometric (ln(mass) linear in n).

- The real SM charged-fermion masses ARE approximately geometric: ln-mass gaps equal
  to within a factor ~1.3-1.9 for leptons, up-, and down-type. Exactly the structure
  this produces.
- The per-generation factor is tunable by the source size (0.15->1.2x, 0.25->2.2x,
  0.4->18x), reaching the lepton ballpark.
- So the SM's "unnaturally tiny" electron Yukawa (~3e-6) is exp(-O(10)) -- the
  exponential of an order-one overlap, not a fine-tuned small number.

Delivered: the SHAPE (geometric, orders of magnitude from O(1) inputs). NOT delivered:
the exact ratios (leptons 207 then 17 -- a *decreasing* ratio, i.e. mild curvature in
log), which needs the detailed level/overlap structure. The flavour problem is not
solved, but it is reframed: "why tiny Yukawas" becomes "internal levels are evenly
spaced," which is what ladders do.

## The self-consistent soliton (the iteration) -- target #1 demonstrated

`self_consistent.py` + `figures/pdf/self_consistent.pdf`. The Hartree loop (fill
levels with Pauli -> density sources the four-fermion sigma field -> that field is
the well -> iterate) CONVERGES (residual ~1e-6 in ~40 iterations) to a genuine
self-bound soliton: bound levels (E < m0), a chiral-restored core (M(0) flips sign),
a finite mass. So Part II target #1 -- the existence the first pass had to ASSUME --
is now demonstrated, not assumed. The gate is open.

S did NOT move to a reliable number: this non-chiral (massive-sigma, Walecka-type)
model gives the soliton but not the symmetry-breaking f_pi and the composite V/A
current correlators that a trustworthy S needs. So: existence moved a lot (the big
one, everything funnels through it); S is still the owed make-or-break, with the
degeneracy result saying only that the *direction* is favourable.

## The chiral soliton: configurational mass, and the limit on S

`chiral_soliton.py` + `figures/pdf/chiral_soliton.pdf`. The symmetry-breaking (linear
sigma / Friedberg-Lee) self-consistent soliton: a real condensate v = f_pi, the
fermion's constituent mass g v, and a self-bound bag.

Delivered (feeds the mass insight):
- a REAL f_pi = v (the order parameter), set by the condensate -- no longer a finger
  on the scale;
- CONFIGURATIONAL MASS: the observable soliton mass is ~94% field + bag energy and
  only ~6% constituent g v -- exactly the proton (~99% binding);
- mass scales linearly with v and vanishes as v -> 0: the condensate CREATES the
  mass. So mass is bound field energy, not a fundamental conserved label; only the
  total energy-momentum (and the charges B-L, Q, spin) is conserved -- which is the
  SAME thing the bounce keeps through the membrane (Part I). Particle-mass insight and
  cosmology close into one statement.

NOT delivered -- S: even with a real f_pi, the S proxy comes out ~10 (unphysical)
because M_V ~ the fermion gap is the WRONG scale. A reliable S needs the composite
vector/axial CURRENT correlators (a further step). So the chiral model cracked the
MASS question, not S. S remains the owed make-or-break; the only directional info is
that degeneracy/walking pushes it favourably.

## The S computation (leading order) -- the make-or-break, honestly

`electroweak_S.py` + `figures/pdf/electroweak_S.pdf`. The constituent-fermion loop
(mass M=gv from the chiral condensate): the V-A spectral difference is positive and
UV-finite (rho_V - rho_A ~ M^2/s), so S > 0. The leading value is the standard
S = N_c/(6pi) per doublet:

- QCD-like (N_c=3, 1 doublet): S = 0.16   (> LEP bound 0.1: graveyard)
- 2 doublets:                  S = 0.32
- Pati-Salam (N_c=4):          S = 0.21   -- extra colour makes it WORSE, as warned.

Escape (S < 0.1) needs M_V/f_pi > ~14 vs QCD's ~8: a ~1.7x 'walk'. That reduction is
the SUBLEADING resonance/meson-loop (walking) correction -- the V and A composites
rearranging so their spectral functions cancel -- which needs the RPA correlators
this calculation does not contain (a lattice-scale problem).

VERDICT: leading-order S from the model is in the graveyard; the escape is genuinely
UNDECIDED. Degeneracy/Pauli says only that the direction is favourable. So Part III
can still die here -- exactly its own stated 'live possibility of being wrong'. The
honest scorecard after the whole soliton thread: existence YES, Regge YES, generations
(shape) YES, configurational mass YES (+ the bounce bridge); S still the coin-flip.

## The Fierz probe (cheap decisive test before the RPA) -- GREEN

`fierz.py` + `figures/pdf/fierz.pdf`. The electroweak S of a composite sector is
killed or saved by the vector/axial resonance splitting M_A/M_V, which is set by the
four-fermion couplings in the meson (exchange) channels, G_V and G_A. So the
make-or-break reduces to one ratio, G_A/G_V, fixed by the interaction's Lorentz
structure.

Fierzed the torsion (Hehl-Datta, axial-axial) term into the exchange channels with
explicit Dirac matrices, VALIDATED exactly against the V-A self-Fierz identity and
Fierz^2 = identity:

    G_S = +1/4,  G_P = -1/4,  G_V = +1/2,  G_A = +1/2,  G_T = 0
    => G_A / G_V = 1  EXACTLY.

The torsion interaction couples the VECTOR and AXIAL meson channels EQUALLY -- forced
by the axial-axial structure, not tuned. That is the structural prerequisite for
walking (M_A -> M_V, where rho_V and rho_A cancel and S -> 0), the OPPOSITE of a
QCD-like sector (G_A != G_V, M_a1 > M_rho, S ~ 0.3, the graveyard). Also G_S = +1/4 > 0
(attractive: drives the chiral condensate the soliton forms) and G_T = 0.

So the cheap probe is GREEN: the framework has a real, structural shot at escaping the
S graveyard that generic technicolor lacks. NECESSARY, not sufficient -- the residual
splitting from the chiral-breaking loops needs the full RPA, which is now clearly
worth building. This flips the leading-order (constituent-loop) verdict's pessimism:
the leading S is the graveyard, but the torsion-specific resonance structure points
hard the other way.

## The RPA (torsion vs QCD-like) -- relative advantage confirmed, absolute owed

`rpa.py` + `figures/pdf/rpa.pdf`. RPA-dressed S = (Nc/6pi) x [resonance enhancement of
the V-A correlator], anchored to the leading value at G=0; walking = large Lambda/M
(small chiral breaking, where Pi_A -> Pi_V).

Result: in the walking limit the loops equalise and the torsion's EQUAL couplings
(G_A=G_V, from the Fierz) keep S bounded while a QCD-like sector (G_A < G_V) blows up
-- S_torsion/S_QCD falls from ~1 (strong breaking) to ~0.36 (deep walking). So the
Fierz direction is confirmed by an independent calculation: the torsion structure is
on the right side of S, increasingly so as it walks.

HONEST LIMIT: this RPA is ONE-SIDED -- it has the vector-resonance enhancement (which
RAISES S, the technicolor problem) but not the full axial/Weinberg-sum-rule catch-up
(which LOWERS S). So the ABSOLUTE values are overestimates and never drop below the
leading ~0.16. It establishes the RELATIVE advantage (torsion << QCD-like), NOT the
absolute escape S<0.1. A trustworthy absolute S needs the full chiral RPA (proper sum
rules, the a1 catch-up) -- likely lattice. So after Fierz + RPA: the torsion is
FAVOURABLY PLACED on S (two independent calcs agree), but the absolute S<0.1 is still
genuinely owed and undecided. Part III is more hopeful than the leading-order
graveyard suggested, not yet safe.
