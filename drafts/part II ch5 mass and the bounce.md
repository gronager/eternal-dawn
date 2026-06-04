# Part II, Chapter II.5 — Mass is Configuration, and What the Bounce Keeps

*Draft chapter for Part II (the Solitonic Origin of Matter), kept separate from the
main Part II file. Register and structure review. This is the capstone of Part II:
it takes the soliton picture's account of mass to its conclusion and finds that the
conclusion is a statement Part I already made from the other end. The two halves of
the book meet here.*

-----

**Status.** The central claim of this chapter — that rest mass is bound field energy,
created by a condensate, and therefore not separately conserved — is **forced** by
the soliton picture and is now **computed**, not asserted: a self-consistent chiral
soliton (`sims/.../chiral_soliton.py`) returns a mass that is ∼94% field-and-bag
energy and only ∼6% constituent, and that scales linearly with the condensate,
vanishing when the condensate is switched off. The further claim — that this is the
*same* conservation law the bounce of Part I obeys — is **structural and exact at the
level of which quantities are conserved**, and is the chapter's reason to exist. What
remains **owed** is the precise quantitative match across the membrane (the subject of
Part IV's bounce simulations). Nothing here is speculative in the way Part III is;
this is the soliton picture's most secure consequence, and it happens to be its most
far-reaching.

-----

## The question mass actually poses

Chapter II.3 said most of a composite's mass is the energy of its constituents'
confinement, not the weight of the constituents — the proton's ∼938 MeV against its
quarks' sub-ten. That was offered as a feature the program inherits for free. This
chapter takes it seriously as a *definition* and asks what follows. If the soliton's
mass is the energy of a bound field configuration, then "mass" is not a substance a
particle is made of. It is a number that describes a *standing pattern* of the
gravity–torsion field. And a number that describes a pattern is a very different kind
of thing from a number that is conserved.

The computation makes the point concrete rather than rhetorical. Solve the chiral
(symmetry-breaking) soliton self-consistently — the fermions fill the levels, their
density melts the condensate in a core, the condensate is the well that binds them,
iterate to a fixed point — and decompose the converged mass. It comes out, for a
representative bag, as roughly six percent fermion level-energy, and the rest field
gradient and bag (vacuum) energy: ∼94% of the observable mass is *field*, not
constituent. This is not an analogy to the proton; it is the same accounting,
reproduced from the soliton's own equations. The "particle" weighs almost nothing of
itself. What weighs is the configuration it holds the field in.

## Mass is switched on by the condensate

There is a sharper way to see that mass is configurational, and the chiral model
shows it directly. The condensate — the vacuum value $v$ of the $\sigma$ field, which
is the order parameter $f_\pi$ — is what gives the fermion any mass at all: the
constituent mass is $g v$, and the soliton's whole mass scales *linearly* with $v$.
Turn the condensate down toward zero and the mass goes with it, to zero. There is no
residual, no irreducible lump of "rest mass" left behind when the condensate is off.
Mass is not a property the fermion carries into the vacuum; it is something the vacuum
*does to* the fermion, and the amount it does is set entirely by how strongly the
vacuum has condensed.

This also disposes of the worry that the mass scale is a finger on the scale. It is
not put in by hand. A condensate scale is generated *dynamically*, by dimensional
transmutation: a coupling run down from a high cutoff reaches strong coupling and
condenses at a scale exponentially below the cutoff — the same mechanism by which QCD
makes its ∼200 MeV from a theory with no small numbers in it. So the one scale that
sets every mass is itself derived, exponentially, from the gravity–torsion coupling,
and the smallness of particle masses against the Planck scale is the exponential of an
order-one number rather than a tuning. The scale is earned, not inserted.

## The consequence: mass is not conserved; the combined quantity is

Here is the step that matters. If mass is the energy of a standing field
configuration, and configurations can be assembled and dissolved, then **mass is not a
separately conserved quantity.** What is conserved is the total energy–momentum
$P^\mu$ of the gravity–torsion field, together with the genuine charges — electric
charge, the difference $B-L$, angular momentum. Mass is the part of $P^\mu$ that
happens to be standing still in a knot, and standing-still-ness is not a conservation
law.

Ordinary physics already says this, in the four words $E = mc^2$ and in every reaction
where the masses do not add up. Beta decay is the cleanest case, and it reads
transparently in the soliton picture: a neutron at 939.6 MeV becomes a proton at
938.3, an electron at 0.5, an antineutrino, and a spray of kinetic energy — and the
rest masses on the two sides simply do not balance. Only $P^\mu$ balances, and only
the charges ($B-L$, electric charge) are carried across unchanged. In the soliton
reading this is one quark-soliton rearranging ($d \to u$) while the energy released
materialises a lepton soliton–antisoliton pair; the *masses* are reshuffled freely
because they were only ever configurations, while the conserved labels ride through
untouched. The neutron's mass was never a thing that had to be conserved; it was the
price of holding three solitons in a particular pattern, and the reaction simply
re-priced them.

So the honest statement of what mass is, is also a statement of what it is not: it is
not a conserved Casimir stamped on a particle for all time. It is the standing energy
of a field pattern, freely convertible to the field's unbound (kinetic, radiative)
energy, with only the total and the charges held fixed.

## Where Part I already said this

And now the chapter's reason for existing. Run the same logic from the cosmological
end. Part I's bounce dissolves everything that falls into it — "not a torn teddy bear,
not a rusty car, not an atom, not even an atomic nucleus" — melting all of it into a
featureless soup at the bounce density, and passing across the membrane *only* the
conserved charges: energy–momentum, $B-L$ (sphaleron-locked), angular momentum,
electric charge. Structure does not cross. The child is rebuilt from elementary field
content on the far side.

Set the two statements side by side:

> *Part II:* mass is configuration; it is not conserved; only $P^\mu$ and the charges
> are.
>
> *Part I:* the bounce dissolves configuration; it conserves only $P^\mu$ and the
> charges.

They are the same statement. The quantity that is *not* conserved in beta decay is the
quantity that is *dissolved* at the bounce — configurational field energy, i.e. mass.
The quantities that *are* conserved in beta decay are exactly the quantities that
*cross* the membrane. The membrane is not doing something exotic to mass that ordinary
reactions do not; it is doing the *same* thing — re-pricing configuration while holding
the invariants — only total, and at the scale of a universe. A beta decay is a bounce
in miniature: a small re-configuration that keeps $P^\mu$ and the charges and discards
the rest. The bounce is beta decay at cosmological scale: the maximal re-configuration,
keeping the same short list and discarding everything else.

This is why the dark sector and the matter asymmetry are the things that survive a
bounce (Part I) and the things that are protected charges (Part II), while the
particular *masses* of the parent's matter mean nothing to the child: masses are
configuration, and configuration does not cross. The framework did not need two
separate accounts of what is fundamental — one for particles, one for cosmology. It
needs one, and the one it needs is the same at both ends: **the invariants are
$P^\mu$, $B-L$, electric charge, and spin; everything else, mass included, is
configuration the field adopts and abandons.**

## What is forced, what is owed

- **Forced (from the soliton picture, and computed):** that the soliton's mass is
  overwhelmingly bound field energy, not constituent; that it is created by the
  condensate and vanishes without it; and therefore that mass is not a separately
  conserved quantity, only $P^\mu$ and the charges are.
- **Structural and exact (the bridge):** that this conserved short list is identical
  to the list the bounce of Part I carries across the membrane — so the
  not-conserved-ness of mass in a reaction and the dissolution of mass at the bounce
  are one fact, seen at two scales.
- **Owed (Part IV):** the quantitative match — that a bounce simulation, tracking
  $P^\mu$, $B-L$, charge, and spin across the membrane while configuration is
  dissolved, reproduces the inheritance Part I asserts; and, on the particle side,
  the precise mass spectrum from the overlaps (Chapter II.4) whose *scale* is the
  transmuted condensate of this chapter.

In one breath, in the book's idiom: *mass is not a substance the universe is made of
and must conserve; it is a configuration the gravity–torsion field adopts, switched on
by a condensate and freely re-priced in every reaction — and the bounce is simply the
most complete re-pricing there is, keeping the same four invariants a beta decay keeps
and dissolving, as it must, everything that was only ever a pattern.* That sentence is
the spine that joins Part I to Part II, and it is the truest thing the soliton picture
has to say.
