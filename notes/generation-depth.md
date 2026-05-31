# Where does our universe sit in the supraverse tree?

Done properly (was postulated/summarized in Ch.5; now derived as a bounded
branching process). Code: `cartasis_sims.population`; figure
`figures/scripts/ch05_population.py` -> `figures/pdf/generation_depth.pdf`.

## The model

Birth: OGUs (== BHU_0) nucleate from the void at a steady rate (infinite,
stationary void). Growth/branching: each viable universe spawns N viable
children -- holes that exceed internal M_crit AND clear the chiral-vortical birth
filter |omega/T| > eta_min/C. A fraction p stay fertility-preserving, so the
effective per-generation reproduction is m = N*p. Universe count at depth n ~ m^n
up to a viability truncation D:

    P(n) = m^n / sum_{k=0..D} m^k.

Three regimes:
- m < 1 (subcritical): lineage dies out, most universes shallow.
- m = 1: scale-free.
- m > 1 (supercritical): population piles at the deepest viable generation; a
  typical universe is DEEP.

Which regime holds is Einstein-Cartan microphysics (N and p). Not yet pinned.

## The robust anchor (independent of m): we are NOT an OGU

We observe dark matter AND dark energy. In this framework both need a parent:
- DM = parent material projected through the membrane;
- DE = parent's ongoing accretion.
An OGU has no parent => no DM, no DE. We have both => **n >= 1.**

This is the one clean conclusion. `population.we_are_ogu` / `min_generation_*`
encode it; the posterior zeroes the n=0 bin.

## Where it connects to data (the "range of Cs" + DM/DE)

- C (chiral-vortical coupling): each universe's |eta| ~ C|omega/T| set per bounce;
  the viable-bounce vorticity distribution + our measured eta (baryon-to-photon
  ~6e-10) constrain C and the birth-filter threshold -> feed N (the branching
  ratio) -> the depth prior.
- DM/baryon ratio (~5.4) and DE w(a) (DESI) constrain the PARENT's mass and
  accretion state -> additional handle on our depth.

So: population distribution gives the prior over n and a range of C; DM&DE give
the floor (n>=1) and will narrow the parent properties. Expect to update the
wallpaper (Ch.10) as these tighten -- drawn depth there is illustrative.

## Honest status

- DERIVED: the branching structure of P(n); the n>=1 floor (robust).
- ILLUSTRATIVE: the specific m, D, structure_decay in the figure; the
  super-critical pile-up sits at the hard cap D (a soft viability cutoff rounds
  it).
- TODO: get N from C + the bounce-vorticity distribution; get the viability decay
  p from per-generation mass/structure budgets; turn DM/baryon + w(a) into a
  quantitative depth constraint.
