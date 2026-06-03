# Part II — The parameter count: banking on S buys most of the Standard Model

*Draft chapter for Part II. Markdown, working register. This is the bookkeeping behind the
claim that tying the forces and masses to the gravity–torsion field puts Eternal Dawn at
**far fewer parameters than the Standard Model**. It is written conditionally — "if we bank
on S" (the electroweak make-or-break) and the lattice targets deliver — and it keeps an
honest column for what is owed. Reproducible from `parameter_count.py` and
`figures/pdf/parameter_count.pdf`.*

-----

## The Standard Model's bill

The Standard Model is extraordinarily accurate and extraordinarily *expensive in inputs*.
Counting the free parameters that must be measured and put in by hand:

| block | count | what it is |
|---|---:|---|
| gauge couplings $g_1,g_2,g_3$ | 3 | the three force strengths |
| Higgs sector $v,\,m_H$ | 2 | the electroweak scale and the Higgs mass |
| charged-fermion masses | 9 | 6 quark + 3 lepton masses (the Yukawa sprinkling) |
| CKM mixing | 4 | 3 angles + 1 CP phase |
| $\theta_{QCD}$ | 1 | the strong-CP angle |
| neutrino masses + PMNS | 7 | 3 masses + 3 angles + 1 phase |

That is **19** without neutrino masses, **26** with them — and the great majority (the 9
fermion masses, the 4 CKM mixings, the 7 neutrino parameters: 20 of the 26) are *pure
insertion*. Nothing in the SM explains why the electron weighs what it does; the number is
measured and typed in. The Yukawa sector is, structurally, a list of empirical constants.

-----

## Eternal Dawn's bill — if we bank on S

The ED particle sector inherits exactly the base constants the rest of the book already
uses: $G,\hbar,c$. It introduces **no new dimensionful input of its own**, because the one
scale the particle sector needs — the condensate $v=f_\pi$ — is *dynamically generated* by
dimensional transmutation of the torsion coupling (just as $\Lambda_{QCD}$ is generated in
QCD), not inserted. So in the ideal, fully-delivered case the particle sector adds **zero
new fundamental free parameters**, and the SM's blocks become *derived* quantities:

- **Gauge couplings (3) — derived.** The strong sector's colour structure is *forced* by
  the three-valued Pauli label (the colour factors equal QCD's, rigorously — `fierz.py`,
  `color_force.py`). The weak/electroweak structure comes from the *one* torsion
  four-fermion coupling, whose Fierz content gives $G_A=G_V$ (walking-favourable). The
  strengths trace to $G$ through the connection. Three couplings $\to$ one field's
  properties.
- **Higgs sector (2) — derived.** There is no elementary Higgs; the condensate is
  composite, and its scale $v=f_\pi$ is generated, not free.
- **9 fermion masses — derived, numbers owed.** Masses are wavefunction *overlap*
  integrals (Yukawa-as-overlap); the geometric generation ladder is shown
  (`generations*.py`), the *exact ratios* are owed to the lattice (punch-list L4).
- **4 CKM mixings — derived, owed.** Mixings are overlaps between soliton generations.
- **7 neutrino parameters — derived, owed.** Same overlap mechanism, plus the CPT-mirror
  structure of Part I.
- **$\theta_{QCD}$ (1) — owed.** Strong CP is *not yet addressed* — an honest residual.

So the structural claim is stark: **~19–26 inserted parameters $\to$ 0 new fundamental
ones.** Where the SM has a list of empirical masses, ED has a single field whose phase
diagram (the genesis chapter) is supposed to *produce* that list. That is the sense in which
banking on S — taking the electroweak sector to be the strongly-coupled gravity–torsion
sector — collapses the parameter count.

-----

## The honest residuals (the gap between promise and proof)

This is a chapter about *structure*, and structure is not proof. The load-bearing caveats:

1. **The numbers are owed to the lattice.** The count above says every SM parameter *can be*
   a consequence; it does *not* yet produce a single one to compared precision. $S<0.1$
   (the electroweak make-or-break, punch-list L3), the exact mass ratios (L4), and even
   *whether the vacuum confines* (L1–L2) are all non-perturbative and unsettled. The
   parameter reduction is real **in structure**, a promise the structure supports — not a
   delivered result.
2. **The gauge group / number of colours is a choice, not yet a derivation.** Three colours
   is taken as the Pauli minimum (it reproduces the observed channel structure and gives the
   lower $S$); four (Pati–Salam) is an elegant but unforced, $S$-worsening option.
3. **$\theta_{QCD}$ is simply not addressed.** Strong CP remains an open residual, listed
   rather than hidden.

The figure (`parameter_count.pdf`) shows both ledgers side by side: the SM's stacked 26,
every block inserted; the ED column, every block traced to $G,\hbar,c+v$, colour-coded by
status (derived / forced / owed). The honest headline is the one in the figure's title:
**far fewer parameters than the SM — the structure supports zero new fundamental inputs;
the numbers are owed to the lattice.**
