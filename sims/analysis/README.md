# Lattice analysis: live health + early spectroscopy

Two standalone scripts (numpy + matplotlib only) to (a) confirm the L4 HMC is sound
and (b) pull m_pi, m_N, m_N/m_pi off the configs you already have -- you do **not**
need to wait for trajectory 1500. Plots land in `sims/output/lattice/`.

## 1. HMC health (run on the live log)

```bash
python sims/analysis/hmc_health.py path/to/RUN.log
# pin the thermalisation cut yourself instead of auto:
python sims/analysis/hmc_health.py path/to/RUN.log --therm 250
# a whitespace table instead of a Grid log (columns: traj dH acc plaq):
python sims/analysis/hmc_health.py ensemble.txt --table
```

Checks: acceptance (want 0.70-0.90), `<exp(-dH)> == 1` (the unbiasedness test),
`<dH>` and the dH histogram, plaquette plateau + first-half/second-half drift test,
and the plaquette integrated autocorrelation time (tau_int -> traj per independent
config). Exit code 1 if any check WARNs. If it parses 0 records, the log format
differs from stock Grid -- paste ~20 lines and the regexes at the top get a one-line fix.

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
