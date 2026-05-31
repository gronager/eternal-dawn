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

## OGU size distribution (birth vs growth) -- DERIVED

Births favour small (P ~ exp(-M/M0)); runaway growth Mdot ~ M^g transports mass
up. Stationary continuity in mass space => upward flux = birth rate => 

    n(M) ~ 1/Mdot ~ M^{-g}   (power law; g=1 Eddington, g=2 Bondi),

cut at high mass by depletion. Many small, few large -- derived, not assumed.
A low-mass viability cut M_vis (enough for astrophysics + BHs) selects observer
universes, which pile up just above M_vis. => we sit near the viability edge, NOT
the rare massive tail. This IS the old "near the optimum" conjecture, now derived.
(cartasis_sims.population.ogu_mass_density / observer_mass_density; fig left panel.)

## Spin distribution (birth vs viability) -- DERIVED shape, amplitude open

Seed vorticity omega sets inherited asymmetry eta ~ C|omega|/T. Sign -> chirality
(even => half matter/half antimatter). Magnitude -> purity. Random fluctuation:
P_birth(omega) Gaussian at 0 (low spin likeliest). BUT |eta|<eta_min => annihilates
=> sterile, evaporates (hellish, no descendants). Viability: |omega|>omega_min.

  Opposite to intuition: no spin = most likely but most hellish; high spin = rare
  but pure. Observer distribution P_obs ~ P_birth * productivity * 1[|omega|>omega_min]
  peaks JUST ABOVE threshold: viable universes are low-spin, low-purity, barely
  past sterile.

Confront with us: our eta ~ 1e-9..1e-8 (small) puts us right at that low-spin edge,
where most observers should be. Predicts small but NONZERO net rotation -> compare
to Shamir JWST galaxy-spin handedness + CMB axis (Ch.8). Amplitude needs C, sigma
(open); shape robust. (population.spin_birth_pdf / spin_observer_pdf / sterile_fraction.)

## Chirality is INHERITED (reconciliation, this turn)

Big change: resolved an existing Ch4/Ch5 contradiction in favour of INHERITANCE.
- Selection happens ONCE at the OG seed (no parent to inherit from; vorticity sign
  picks matter vs antimatter, CPT-symmetric => half/half OGUs).
- Inherited down each lineage: a matter OGU has ONLY matter descendants, antimatter
  OGU only antimatter. Lineages chirally PURE, no mixing. Purity ~1:1e8 set at top,
  amplified (never flipped) by mean-field at later bounces.
- "Hell upward" resolved by inheritance (no progressive symmetrization; climb ends
  at the OG seed), NOT by per-bounce genesis. Ch4 "amplifier-not-source" was already
  saying this; Ch5's old "fresh per bounce, random sign" text was the contradiction
  and is now fixed. Wallpaper updated to match (half/half OGUs, pure lineages).

## Honest status

- DERIVED: the branching structure of P(n); the n>=1 floor (robust).
- ILLUSTRATIVE: the specific m, D, structure_decay in the figure; the
  super-critical pile-up sits at the hard cap D (a soft viability cutoff rounds
  it).
- TODO: get N from C + the bounce-vorticity distribution; get the viability decay
  p from per-generation mass/structure budgets; turn DM/baryon + w(a) into a
  quantitative depth constraint.
