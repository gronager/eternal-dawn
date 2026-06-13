# Lattice analysis: live health + early spectroscopy

Two standalone scripts (numpy + matplotlib only) to (a) confirm the L4 HMC is sound
and (b) pull m_pi, m_N, m_N/m_pi off the configs you already have -- you do **not**
need to wait for trajectory 1500. Plots land in `sims/output/lattice/`.

## 1. HMC health (run on the live Grid log)

```bash
python sims/analysis/hmc_health.py lattice/out/dyn_L16x64_m-0.75/stream1/hmc.log
# pin the measurement thermalisation cut (in trajectories) if auto looks off:
python sims/analysis/hmc_health.py .../hmc.log --therm 100
# override beta/volume if the banner isn't in your slice of the log:
python sims/analysis/hmc_health.py .../hmc.log --beta 5.6 --volume 262144
```

Reads Grid's `GridLogHMC` summary format directly:
`Total H after trajectory ... dH = X`, the `Skipping Metropolis test` thermalisation
phase, and `Metropolis_test -- ACCEPTED/REJECTED`. Plaquette is reconstructed from the
Wilson gauge action (`S [1][0] H`) as `1 - Sg/(beta*6*V)` -- beta and volume are
auto-read from the banner (`beta=...`, `Full Dimensions : [...]`).

Reports: thermalisation length (Metropolis-skip count) and where the plaquette plateau
sets in; acceptance and `<exp(-dH)>==1` over the PRODUCTION trajectories only; `<dH>`;
plaquette plateau mean + first/second-half drift; and `tau_int(plaq)` -> how many traj
between independent configs (so you know the right config-save stride and your N_eff).
Exit 0 healthy, 1 if a WARN fires, 2 if the parse is too thin to judge.

## 2. Spectroscopy plateau (run on thermalised configs, now)

```bash
# auto-pick the plateau windows first, eyeball the effective-mass tables/plots:
python sims/analysis/spectrum_plateau.py --pion pion.dat --nucleon nucl.dat --nt 64

# then pin the windows for the final number (cosh for the meson, exp for the baryon):
python sims/analysis/spectrum_plateau.py \
    --pion pion.dat --nucleon nucl.dat --nt 64 \
    --pion-type cosh --nucleon-type exp \
    --pion-window 16 28 --nucleon-window 8 18
```

Input is any whitespace/CSV text in one of three layouts (auto-detected):
`cfg t C` per row, or `t C` in repeated config-blocks, or a `configs x Nt` matrix.
Output: `m_pi a`, `m_N a`, and the correlated-jackknife `m_N/m_pi` (the number s_T
keys on; constituent-counting target 3/2). Feed that ratio into the existing `run/10`
ground-state line for s_T itself -- these scripts stop at the lattice observables.

Tip: ratios converge faster than absolute masses, so a few tens of thermalised configs
already bracket `m_N/m_pi`. Start a parallel measurement job on the saved configs while
HMC keeps generating.

## 3. Chiral extrapolation (valence-mass scan)

After running `06` at several valence masses and saving each `baryon_raw.dat`
(e.g. `baryon_raw_m-0.75.dat`, `…m-0.5.dat`, `…m-0.3.dat`):

```bash
python sims/analysis/chiral_extrap.py \
    baryon_raw_m-0.75.dat baryon_raw_m-0.5.dat baryon_raw_m-0.3.dat
# explicit masses / windows:
python sims/analysis/chiral_extrap.py f1 f2 f3 --masses -0.75 -0.5 -0.3 --T 64 \
    --pion-window 8 21 --nucleon-windows 10,20 12,24 16,26
```

Fits m_pi a (cosh) and m_N a (forward log, over several windows -> the window systematic
that dominates the baryon error) per ensemble, then extrapolates: the **GMOR line**
(m_pi^2 vs valence mass -> critical mass m_crit) and the **chiral nucleon**
(m_N vs m_pi^2 -> m_N^0, the chiral-limit torsiton mass). Lattice units; multiply by 1/a
(w0/r0 from `run/03`) for MeV. Masses are inferred from the `m-*.dat` filenames.
