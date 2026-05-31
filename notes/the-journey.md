# The two-sided journey: nesting, dark matter, and time

Working note. Granting universe = black-hole interior (CPT-even bounce,
Chapter 9 + bounce-transform note), fly the Tardis up and down the nesting and
ask what the physics forces. Three threads — depth, dark matter, time — and they
interlock into clean endpoints in both directions. Robust claims and speculative
ones are labelled.

## 0. The Copernican hook

The recursion is the one part of the framework that is *directly observable*: we
cannot see our own Cartasis membrane, but we see other universes' membranes from
outside — they are the black holes in our sky. Their statistics (how many, what
masses, what they do to gravity and clocks) are a trace of the same process that
made us. That is the Copernican lever.

## 1. Zoom IN — how deep does the tree go? (ROBUST: shallow)

Mass is **not** preserved per level. A child universe's total mass-energy is the
mass of the parent black hole it bounced from, and that is a tiny fraction of the
parent universe's mass:

- our universe: M0 ~ 9.2e52 kg;
- largest known black hole (TON 618, ~6.6e10 Msun): ~1.3e41 kg;
- so the biggest possible child/parent mass ratio is f_max ~ 1.4e-12.

Iterating along the *most massive* branch, m_k ~ M0 * f_max^k:

| level | mass (kg) | what it is |
|------:|-----------|------------|
| 0 | 9.2e52 | our universe |
| 1 | 1.3e41 | a one-large-galaxy-worth universe |
| 2 | 1.8e29 | ~0.1 Msun of mass-energy as an entire cosmos |

A level-1 child holds about one galaxy's mass as its *whole* universe — it can
maybe make a few stellar black holes, not supermassive ones. A level-2 child is
sub-stellar in total: it cannot form black holes at all. So along the richest
branch the tree **terminates in ~1–2 generations downward**, and typical
(stellar-mass) branches terminate immediately.

**Result:** the nesting is *not* an infinite self-similar fractal downward. It is
a rapidly terminating tree, because each child is ~12+ orders of magnitude
lighter than its parent and the Standard-Model floor for forming a black hole
(needing roughly a galaxy's mass) is only ~12 orders below us. "When do we stop
seeing new black holes?" — almost immediately, going in.

## 2. Zoom IN — what happens to dark matter? (CRUX: the projection law)

If dark matter is the parent's gravitational influence felt through the membrane
("outside" gravity), then the dark-to-visible ratio of a universe is set by how
much outside mass it feels relative to its own visible mass.

Naive scaling: a child of mass m sits inside a parent of mass M ~ m/f_max ~ 7e11 m.
If it felt the *whole* parent, its dark/visible ratio would be ~7e11 : 1 — almost
purely dark. But **we observe ~5.4 : 1, not ~1e12 : 1.** That single fact is a
strong constraint:

> The membrane projects only a tiny, specific fraction (~1e-11) of the parent's
> mass into the child — or the felt quantity is not the total parent mass at all
> (e.g. only mass within a membrane-scale region, or the parent's own dark
> component). Pinning this **projection law** is the key open calculation; the
> observed modesty of dark matter is its first data point.

Two candidate laws, with opposite recursion fates:

- **(a) Depth-amplifying:** felt-outside grows relative to own mass as you nest
  down → dark fraction rises level by level → deep children are dark-dominated,
  which suppresses their structure formation and *reinforces* the downward
  termination of §1. (Speculative; needs the law.)
- **(b) Fixed point:** if the nesting is self-similar in mass ratio, the
  dark/visible ratio is the same ~5:1 at every level — and we measure ~5:1
  because that is the universal fixed point of the Cartasis recursion. Elegant,
  Copernican, but requires self-similar mass ratios, which §1 breaks at the
  downward floor. (Speculative; appealing.)

Either way one endpoint is **robust**: an OG universe nucleated directly from the
void has *no parent*, hence no outside influence, hence **no dark matter**
(within the framework's own definition of DM). See §4.

## 3. Time — frozen going in, hidden going out (ROBUST: from the OS solution)

The bounce does not reverse time (Janus point), but gravitational time dilation
is real and we already computed it (the OS "frozen star", Chapter 9): the surface
of a forming black hole freezes at the horizon as seen from outside — exterior
Schwarzschild time t -> infinity and redshift 1+z -> infinity as R -> r_s.

Consequences for the journey:

- **Zoom in:** a child universe's *entire internal history* (billions of its
  years) plays out behind its horizon. From our vantage we never watch it tick
  slowly forever — we watch its formation redshift to black and then it is
  causally sealed. Its whole cosmos is "frozen" at formation from outside.
- **Zoom out:** symmetrically, *our* entire 13.8 Gyr is sealed behind our horizon
  from the parent's exterior — our universe is, to the parent, a black hole that
  formed and went dark. The parent's clock is not slowed by us at all.

So the clock relationship across a membrane is the divergent interior↔exterior
map of the OS solution, not a finite rate ratio. **To compute:** continue that
map through the bounce to relate baby cosmic age to parent exterior time
explicitly (we have all the pieces in `os_collapse.py`).

## 4. Zoom OUT — the logical end (ROBUST endpoint: OG = dark-matter-free)

Going up, universes get more massive (~1/f_max ~ 1e12 per level) until the chain
reaches an **OG universe**: one nucleated directly from the primordial void, with
no parent. By the dark-matter mechanism, an OG universe has nothing "outside" to
project gravity in, so it carries **no dark matter**. The upward chain therefore
terminates at dark-matter-free, void-seeded roots.

This hands dark matter a second job: it is a **depth gauge**. A universe's
dark/visible ratio measures how deep in the nesting it sits — zero at the OG
root, nonzero for every descendant. We measure Omega_dm/Omega_b ~ 5.4, so **we
are a descended universe, not an OG root** — and (under law (b)) perhaps at the
recursion's fixed point.

## 5. Does it continue forever, or does dark matter make it diverge?

Neither direction is infinite, and the boundaries are different in kind:

- **Downward:** terminates fast (§1) at the Standard-Model mass floor for black
  hole formation — accelerated if dark matter is depth-amplifying (§2a).
- **Upward:** terminates at OG/void roots (§4), which are dark-matter-free.

Dark matter is the order parameter that runs monotonically between the two ends:
it is the trace of having a parent, vanishing at the top and (plausibly) growing
toward the bottom. The recursion is a **finite tree, bounded above by the void
and below by the mass floor**, with our ~5:1 placing us in the interior.

## 6. To compute / next

1. **Nesting cascade** (cheap): formalise §1 with the observed black-hole mass
   function; predict the depth distribution and "deepest observable child."
2. **Projection law** (the crux): derive how membrane junction conditions project
   parent mass-energy into the child; check it can give ~5:1, and determine
   whether it is depth-amplifying (a) or fixed-point (b).
3. **Interior↔exterior time map** (medium): continue the OS map through the
   bounce; quantify "a baby's full history = frozen at our horizon."
4. **OG = no dark matter** (theory): state as a sharp, in-principle-testable
   prediction; contrast OG vs descended universes.

## 7. Tee-up to Direction 2 (the void)

§4 dead-ends upward at the primordial void. The companion idea: the void is not
structureless — it is the most thermodynamically stable manifold, and
fluctuations do not just happen, they "land somewhere." That shape, and how a
fluctuation condenses on it, is the genesis question — taken up separately.
