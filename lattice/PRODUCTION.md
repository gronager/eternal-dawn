# Production campaign — the resolved, unitary torsiton bag

The pilot (12³×24, β=5.6, sea −0.5, 48 cfg) validated the whole pipeline but hit its **resolution
limit**: the bag is ~1 lattice spacing, and T=24 starves the two-state fit (RESULTS.md, L4). This
campaign grows the lattice to resolve the bag and pin `s_T`, `g_S`, the chiral trend, and the
generation count. No rebuild is needed — every binary takes `--grid` at runtime.

## Geometry and why

Stay at **β=5.6** (same `a`); the bag is resolved not by finer `a` but by a **lighter quark** (the
bag grows as m_π falls — the pilot proxy already showed this), which a bigger box enables.

| change | buys |
|---|---|
| **T = 48** | backward node moves t=12→24 ⇒ a clean sink at large `t_snk` *and* ~16 τ-points for the two-state fit (the biggest 3-pt win) |
| **L = 16** | finite-volume control, room for extended radial excitations (GEVP), and the bigger volume + dynamical determinant push back the exceptional-config wall ⇒ a lighter sea is reachable |
| **lighter sea** | bigger, *resolvable* bag + the **unitary** chiral point (retires partial quenching) + a slightly finer `a` (sea softens confinement) |

Deferred to a real machine: the **continuum limit** (≥2–3 spacings at finer β). Not in this campaign.

Cost vs the pilot: volume ×(48/24)·(16/12)³ ≈ **×4.7**. The lighter-sea masses cost *more* on top —
the fermion force stiffens near the chiral limit (more MD steps, longer thermalisation, more CG).

## Stage 1 — method check (sea −0.5, heavy, cheap)

Same heavy mass as the pilot, just bigger T,L. The only question: **does the bag resolve and the
two-state fit stabilise once T=48 gives it τ-room?** This isolates the T-effect from the mass-effect
before spending on light masses.

```bash
# generate the ensemble (the long step; watch acceptance in stream1/hmc.log)
L=16 T=48 BETA=5.6 MASS=-0.5 NTRAJ=1500 MDSTEPS=28 SAVE=5 NSTREAMS=1 \
  OUT=out/dyn_L16x48_m-0.5 ./run/07_dynamical_torsiton.sh 2>&1 | tee dyn_L16x48_m-0.5.gen.log
```

Watch `out/dyn_L16x48_m-0.5/stream1/hmc.log`: **acceptance ~0.7–0.9** and **dH ~ O(1)** per trajectory
are healthy. If acceptance < 0.7 raise `MDSTEPS` (→32, 40); if it sits > 0.95 you can lower it to save
time. The plaquette should settle within a few hundred trajectories — that is thermalisation.

Then measure (set `THERM` past where the plaquette settles, `STRIDE` past the autocorrelation):

```bash
# 1. scale on THIS ensemble -> read off r0/a (you'll feed it to 09/10 as R0A)
L=16 T=48 OUT=out/dyn_L16x48_m-0.5 THERM=300 STRIDE=10 ./run/01_puregauge_potential.sh
# 2. spectrum (pion + torsiton); confirm a clean m_N plateau
L=16 T=48 MASS=-0.5 OUT=out/dyn_L16x48_m-0.5 THERM=300 STRIDE=10 ./run/06_baryon_spectrum.sh
# 3. one-body bag proxy
L=16 T=48 MASS=-0.5 R0A=<r0/a from step 1> OUT=out/dyn_L16x48_m-0.5 THERM=300 STRIDE=5 \
  ./run/09_bag_profile.sh
# 4. the genuine 3-pt bag -- SINKT defaults to 3T/8=18 (well before the t=24 node), ~16 tau-points
L=16 T=48 MASS=-0.5 R0A=<r0/a> OUT=out/dyn_L16x48_m-0.5 THERM=300 STRIDE=5 \
  ./run/10_condensate_3pt.sh
# 5. generations (nice-to-have)
L=16 T=48 MASS=-0.5 ITERS="0,10,30,60,100" OUT=out/dyn_L16x48_m-0.5 THERM=300 STRIDE=10 \
  ./run/08_torsiton_gevp.sh
```

**Stage-1 success = the run/10 GROUND-STATE line resolves** (R0 ≳ 1.5a, not below the cutoff) and is
stable as you vary `SINKT` (try 16, 18, 20). If it does, the method is sound and we go to light mass.

## Stage 2 — the physics (lighter sea, unitary, chiral) — **at T = 64**

Two lighter sea masses for a chiral extrapolation, **on a longer lattice, T = 64**. This is the key
change from Stage 1: the sink-time scan (RESULTS.md, L16×48) showed the three-body bag `s_T` *drifting*
0.45→0.56 across sink times because the two-state fit starved for τ-room and `t_snk=20` crowded the
backward node at T/2 = 24. **T = 64 moves the node to 32**, giving the fit the room to *pin* `s_T`
rather than bracket it — the decisive lever, even more than the lighter sea. (The Hasenbusch
convention is **validated** — see below — and is geometry-independent, so it is trusted here.)
Targets are by **m_π a**: tune toward m_π a ≈ 0.5 and ≈ 0.35.

```bash
# launch DETACHED so a dropped terminal can't kill it; run/07 auto-resumes if it ever does.
for M in -0.75 -0.85; do
  nohup env L=16 T=64 BETA=5.6 MASS=$M HASEN="-0.3,0.1" NTRAJ=2500 MDSTEPS=20 SAVE=5 NSTREAMS=1 \
    OUT=out/dyn_L16x64_m$M ./run/07_dynamical_torsiton.sh > dyn_L16x64_m$M.gen.log 2>&1 &
done
disown -a    # or, cleaner: run the whole loop inside `tmux new -s prod` and detach with C-b d
```

**MDSTEPS=20, not 40.** The validation run at m=−0.75 showed MDSTEPS=40 gives |dH| ~ 0.03 at
acceptance ≈ 1 — badly over-integrated. Since |dH| ∝ (1/MDSTEPS)², MDSTEPS≈20 lands |dH| ~ 0.1–0.15
at ~0.95 acceptance: a free ~2× on top of the accelerators. For the stiffer m=−0.85 start at ~24–28
and watch acceptance. Then repeat measurements 1–5 per ensemble. Deliverables: `s_T(m_π²)` and
`g_S(m_π²)` extrapolated to the chiral point on a resolved, unitary, T=64 lattice.

### Speeding up the light-sea generation (the GH200's 96 GB)

`generate_dynamical` now ships two standard dynamical-Wilson accelerators, risk-ordered:

- **Even-odd (red-black) preconditioning — always on, free, exact.** A determinant *identity*, so it
  cannot bias the ensemble (worst case: lower acceptance). ~2×. No validation beyond the usual
  acceptance/dH check.
- **Hasenbusch mass preconditioning — opt-in via `HASEN`, off by default.** Splits the determinant into
  a telescoping product of mass-ratio factors with smoother forces → bigger steps / fewer CG iters,
  another ~2–3×. EXACT *only if* Grid's numerator/denominator ratio convention is as assumed; a
  backwards convention silently samples the **wrong** action while acceptance still looks fine.

```bash
# light sea WITH the accelerators (Hasenbusch masses HEAVIER than the sea, ascending):
L=16 T=48 BETA=5.6 MASS=-0.75 HASEN="-0.3,0.1" NTRAJ=2500 MDSTEPS=40 NSTREAMS=4 \
  OUT=out/dyn_L16x48_m-0.75 ./run/07_dynamical_torsiton.sh
```

**Validation — DONE (✓ certified).** The protocol was: (a) it compiles; (b) a short run shows
acceptance ~0.7–0.9 and dH ~ O(1); (c) cross-check the plaquette against the plain ensemble — they
must agree within errors. Result: the log telescope reads correctly
(det(−0.75)/det(−0.3)·det(−0.3)/det(0.1)·det(0.1) = det(−0.75)); dH is symmetric and sub-0.05; and a
**stationarity test** — the exact production telescope (HASEN=−0.3,0.1) launched from a *plain*-
equilibrated m=−0.5 config — held the plaquette at **0.5514**, dead on the 35-config plain reference
**0.55138(1)** (difference 3×10⁻⁵, ~20× below the fluctuation). A backwards convention would have
sampled a grossly different action and drifted off within a trajectory or two; it did not. Since the
convention is mass- and geometry-independent, this certifies HASEN=−0.3,0.1 for the T=64 light-sea
runs. Even-odd alone (no `HASEN`) needs only (a),(b). Raise `NSTREAMS` (each stream ~2–3 GB at T=64)
to run many chains in parallel on the one GPU — N× the configs per wall-clock until bandwidth saturates.

**Expected wall-clock** (3-day plain baseline): even-odd ~×2 → ~36 h; + Hasenbusch (2 masses) ~×2–2.5
→ **~14 h, already under a day**; + deflation a further ~×2 → ~6–7 h. The first two are the safe,
validated-by-design stack; deflation is the experimental last factor.

### The VRAM filler — low-mode deflation (`DEFLATE`, `-DUSE_DEFLATION`)

`generate_dynamical` carries an **experimental** deflation accelerator (the genuine 96 GB exploit):
hold `N` low eigenvectors of the sea operator resident and deflate the CG. It is **behind a compile
flag** so it can never break the validated even-odd+Hasenbusch binary:

```bash
make -C lattice/src CXXFLAGS="$(grid-config --cxxflags) -DUSE_DEFLATION" generate_dynamical
L=16 T=48 MASS=-0.75 HASEN="-0.3,0.1" DEFLATE=200 NSTREAMS=2 ... ./run/07_dynamical_torsiton.sh
```

It is **exact** (the action solve refines to CG tolerance) and **reversible** (the subspace is frozen
for the run) by construction — so it cannot bias the ensemble, only the speed. Two honest caveats:
the Grid eigensolver API is version-specific (the `USE_DEFLATION` block flags the exact calls to
fill in/verify), and the subspace **goes stale** as the gauge field decorrelates — recompute it on
each checkpoint/restart to refresh, until the in-run trajectory-boundary refresh hook is added (owed).
*Owed follow-up:* mixed-precision CG, and the sustained in-run subspace refresh.

## What lands

- a **resolved** ground-state bag `s_T` from the two-state fit (T=48 τ-room) → in/out of the window, decisively;
- the **unitary** chiral extrapolation of `s_T` and `g_S` (no partial quenching);
- a cleaner generation count (bigger box for the radial tower; nodal operators if we build them);
- everything cross-checked against the validated pilot pipeline.
