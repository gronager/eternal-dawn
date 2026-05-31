# The void and genesis, through the filtering lens

Working note. The upward end of the journey (notes/the-journey.md) dead-ends at
OG universes nucleated directly from the primordial void. This note develops what
the void is and how genesis works once we read the bounce as a CPT-even
continuation with a CVE-style, parity-violating *filter* (Chapter 4 rewrite +
cve-filter note) rather than a flip. The headline: genesis needs no "hell" — the
void makes baryon-asymmetric, viable universes at first crossing. Speculative
parts are labelled.

## 1. What the void is

The void is the infinite quantum vacuum substrate of Chapter 1 (Axiom 1) +
infinity. Two structural facts:

- **It is the low-gravitational-entropy phase.** Smooth, near-zero Weyl
  curvature, no horizons (Penrose's gravitational-entropy ordering). By the same
  Second-Law argument that drives condensation (Chapter 5), the void is the
  *dilute* phase and the foam of universes is the *condensed* phase. The void is
  locally stable (it is the vacuum) but globally metastable toward condensation.
- **It is CPT-symmetric.** No preferred charge, parity, or time orientation;
  averaged over the substrate, equal capacity for matter- and antimatter-biased
  outcomes.

"Most thermodynamically stable manifold" is then sharpened: the void is the
minimum-gravitational-entropy background, and genesis is the Second Law acting on
it. (Whether the void's geometry can be *derived* by extremising some functional
is the deep open question of §6.)

## 2. Genesis: a fluctuation that lands

"Fluctuations don't just happen, they land somewhere." Concretely: a vacuum
fluctuation is a genuine event only when it reaches the Cartasis density and
bounces — that is the OG nucleation. The nucleation rate (Chapter 5) is

    rate(M) ~ exp(-M c^2 tau / hbar) x productivity(M),

so the void preferentially nucleates *small* OG universes (the exponential
punishes large M), balanced by the requirement that the universe be productive
enough to matter. The void is therefore not a sea of equal universes but a
distribution sharply weighted to the smallest viable bounces.

## 3. Genesis through the filter: matter at first crossing

This is the new part. An OG universe has *no parent*, so its asymmetry cannot be
inherited (Mechanism C needs a seed; there is none upstream of the void). It must
be **generated at the OG bounce itself.** And we now have a mechanism that does
exactly that:

- A vacuum fluctuation large enough to bounce generically carries angular
  momentum (vorticity); a non-rotating fluctuation is measure-zero.
- The rotating OG bounce drives chiral-vortical transport (cve-filter note):
  a parity-odd, spin-axis-aligned chirality flux that, with the sphalerons active
  at the ~10^7 GeV bounce, freezes a baryon asymmetry |eta| ~ C |omega/T| into
  the new universe **at birth**.
- The **sign** of the asymmetry is set by the sign of the fluctuation's
  vorticity (which way it spins). Vorticity sign is random over the void, so OG
  universes are born in equal numbers matter- and antimatter-biased:
  **global CPT is preserved, local CPT is broken per universe.**

So genesis = the void condenses rotating fluctuations into bounces, and each
bounce's rotation writes a baryon asymmetry into the child. Matter is made at the
membrane, once, at the very first crossing.

## 4. No hell (the dissolution of the hell-upward problem)

The old worry (Chapter 5): climbing the lineage upward, filter compounding makes
ancestors progressively more symmetric -> near-symmetric, sterile, "hell." The
new picture removes it on two counts:

1. **No flip, no alternation to climb.** With the CPT-even bounce there is no
   matter/antimatter sign-flip per generation (Chapter 4), so there is no chain
   of alternating, progressively-symmetric ancestors to worry about.
2. **Asymmetry is generated per bounce, not diluted along a lineage.** Because the
   CVE filter *creates* |eta| ~ C|omega/T| at *every* bounce (including the OG
   one), a universe's asymmetry is set by its own bounce rotation, not by how far
   it sits from a symmetric ancestor. The upward end is therefore not brimstone:
   it is the void, and the OG universes born from it are immediately
   baryon-asymmetric and viable provided their bounce rotated enough (§5).

This **closes the Chapter 4 open item** "re-derive the hell-upward resolution
without the flip." Resolution: there is no hell because matter genesis is local
to each bounce; the void is symmetric, each child breaks the symmetry itself, and
the population stays globally CPT-symmetric.

## 5. The birth filter: a viability fraction

Genesis has a built-in selection. A universe is viable (leaves matter after
annihilation, forms structure, makes black holes -> has descendants) only if its
birth asymmetry exceeds a threshold, |eta| > eta_min. Since |eta| ~ C|omega/T|,
this is a cut on the bounce vorticity:

    viable  <=>  |omega/T| > eta_min / C.

The void's vorticity distribution (unknown, but plausibly peaked near zero with a
tail) then sets the *fraction* of OG bounces that are viable. Slowly-rotating OG
fluctuations bounce but stay near-symmetric and sterile (they annihilate to
radiation and make no descendants); fast-rotating ones are viable and seed
lineages. This replaces the "hell-upward climb" with a one-step birth filter, and
it is **computable** given a vorticity distribution and eta_min (see §7).

Note the elegance the 41x shortfall buys us: we do NOT need every bounce to be
maximally asymmetric. Even modestly-rotating OG universes can be viable, and the
spin axis they are born with is the same axis that later shows up in their CMB
and galaxy-spin statistics. One axis, from genesis onward.

## 6. What the void does NOT explain (honesty)

Relocating genesis to the void sharpens but does not close the deepest questions:

- **The measure problem.** Counting observers across an infinite void to make
  "we are typical" meaningful is unsolved (shared with all multiverse pictures).
- **The cosmological constant problem.** Why the void is nearly flat rather than
  de Sitter at the QFT cutoff is untouched, arguably made sharper by an infinite
  active vacuum.
- **What fixes the laws.** The Standard Model + Einstein--Cartan are inputs; the
  void does not derive them. The regress is relocated, not terminated.

These belong in the open-questions appendix; the void is genesis *within* given
laws, not an account of the laws.

## 7. To compute / write

1. **Viability fraction** (cheap, illustrative): assume a vorticity distribution
   p(omega/T) and a threshold eta_min; compute the fraction of OG bounces that are
   viable and confirm the population is CPT-symmetric (equal matter/antimatter).
   A small Monte Carlo + figure; clearly parameter-dependent.
2. **OG size distribution peak** (illustrative): combine the exp(-M c^2 tau/hbar)
   suppression with a productivity model; locate the peak; compare to our mass.
3. **Propagate to the monograph:** fold "no hell / genesis through filtering" into
   Chapters 4 and 5 (replace the flip-dominated lineage text), and add the void
   limitations to the open-questions appendix.
