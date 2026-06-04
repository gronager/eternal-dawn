# Part IV — Addendum: the lattice punch-list

*Draft addendum for Part IV (Observations, Simulations, Experiments). Markdown, in the
same register. This is the honest hand-off: the calculations the matter-and-force
programme of Parts II–III has reduced to well-posed non-perturbative problems, specified
precisely enough that a lattice gauge-theory group could pick them up. It is written so
that the boundary between what laptop-scale computation settled and what genuinely needs
a supercomputer is unmistakable.*

-----

## Why there is a punch-list at all

Parts II and III were computed — but with analytic, few-body, and mean-field methods
(nonlinear Dirac solitons, the Hartree loop, Fierz algebra, RPA moments, SU(3) colour
factors, a Nielsen–Olesen vortex). Those reach a remarkable amount: the soliton exists,
mass is configurational, colour structure is forced and matches QCD, generations are a
geometric ladder with a finite cap, and confinement has a candidate mechanism with a
tension set by the condensate. What they do **not** reach is the genuinely
non-perturbative core of a strongly-coupled gauge theory: whether the vacuum confines,
the absolute electroweak S of a near-conformal sector, and the exact mass spectrum. Those
are lattice problems — and lattice gauge theory is one of the most computationally
demanding endeavours in science (Monte-Carlo sampling of the gauge-field path integral on
4-D lattices; leadership-class machines, large allocations, expert teams). So the residue
is not "run it on a bigger laptop"; it is a real research programme. The point of this
addendum is to state it as a *finite, ordered* programme, with success criteria and a
sense of cost, so that the framework's quantitative claims are *consequences to compute*,
not parameters to insert.

A note on the object to be computed. The inter-source potential decomposes as a
short-distance **overlap** (the gravity–torsion field shared in the solitons' tails;
perturbative; screens) plus a long-distance **tube** (the field organised by the
substrate condensate into a flux string; non-perturbative; ~L) — one field with two
limits, the Cornell form V(r) ≈ −κ/r + σL. Lattice is needed for the *tube* (and for the
strongly-coupled spectrum); the *overlap* is already laptop-reachable.

-----

## The targets, in dependency order

**L1 — Does the gravity–torsion vacuum confine? (the cleanest, do first.)**
Put the Part I non-abelian connection on a Euclidean lattice and measure the **Wilson-loop
area law** (and/or a dual order parameter for monopole condensation). *Question:* does the
static potential rise linearly, V(L)=σL, or screen to a constant?
*Success:* a clear area law → confinement, and a string tension σ; *failure:* a perimeter
law / screening → the soliton-confinement picture dies here (a clean falsification).
*Cost:* the *qualitative* confinement question (area vs perimeter) is the **cheapest**
lattice target — modest lattices, no chiral fermions required for the pure-gauge string —
and is the natural first step. It directly tests the `flux_tube.py` assumption (dual
superconductivity) that the laptop computation had to make.

**L2 — The string tension in physical units.** Given L1, extract σ and check the
prediction σ ∝ f_π² (the condensate scale that also sets the configurational mass). *Success:*
a σ consistent with the condensate scale the chiral soliton fixes — confirming "one
substrate, mass and confinement both." *Cost:* moderate (follows L1).

**L3 — The electroweak S of the walking soliton sector (the make-or-break).**
Compute the vector and axial-vector current correlators of the strongly-coupled,
near-conformal sector (couplings G_A=G_V from the Fierz) and extract S. *Question:* does
the walking pull S below the LEP bound (~0.1)? *Success:* S < 0.1 → Part III lives;
*failure:* S in the technicolor graveyard → retreat to Part II intact. *Cost:* **frontier.**
Near-conformal ("walking") dynamics are notoriously hard on the lattice (large lattices,
critical slowing-down); this is the hardest target and the field has not settled the
analogous technicolor question in forty years. The laptop work (Fierz green, RPA relative
advantage) says the framework is *favourably placed*, not safe; only this decides it.

**L4 — The fermion mass spectrum as overlap integrals (the headline prize).**
With lattice soliton wavefunctions in hand, compute the Yukawa couplings as overlap
integrals — one fermion mass ratio, parameter-free, rather than inserted. *Success:* one
ratio matching observation with no tuned input would be the result that makes the whole
programme worth more than its coherence. *Cost:* high (needs the soliton spectrum and the
overlaps at lattice precision). The generation *shape* (geometric ladder) is already
laptop-computed; the exact *ratios* (e.g. leptons 207, 17) are this target.

**L5 — The particle-scale binding density vs cosmological ρ_C (reconciliation with Part I).**
Determine whether the density at which a particle-scale soliton stabilises equals the
cosmological bounce density ρ_C ∼ 10⁵⁰ kg m⁻³, or differs (the bare Einstein–Cartan
coupling naively puts the soliton near the Planck mass unless the relevant density
differs). *Success:* a derived particle-scale density and an explicit relation to ρ_C —
closing the open caveat that the same potential acts in two regimes. *Cost:* tied to L4.

-----

## What is laptop-settled vs lattice-owed (so the boundary is unmistakable)

| Result | Settled (laptop) | Owed (lattice) |
|---|---|---|
| Soliton exists, self-bound | ✅ (`self_consistent.py`) | — |
| Mass is configurational; f_π real | ✅ (`chiral_soliton.py`) | exact masses (L4) |
| Colour charge structure = QCD | ✅ forced (`fierz.py`, `color_force.py`) | — |
| Generations: geometric ladder, finite cap | ✅ shape (`generations*.py`) | exact ratios (L4) |
| Confinement | mechanism + σ∝f_π² *conditional* (`flux_tube.py`) | **does the vacuum confine? (L1–L2)** |
| Electroweak S | direction green (Fierz+RPA) | **absolute S<0.1? (L3)** |

-----

## The dream, and the honest order

The dream is every particle mass to twenty decimals from no free parameters. The honest
order is the reverse of glamour: **L1 first** — the cheapest, most decisive, most
falsifiable (does the vacuum confine, yes or no) — because it is the one that can *kill*
the picture quickly and cheaply, and a framework earns trust by exposing its cheapest
death first. Only if L1 confirms confinement do L3 (the S make-or-break) and L4 (the
masses) become worth their large cost. Every entry above is reproducible from the laptop
code as far as laptop computation reaches; this addendum marks exactly where that reach
ends and the supercomputers begin.
