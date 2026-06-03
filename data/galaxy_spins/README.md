# Galaxy spin catalogue

## gz1_clean.csv  (REAL DATA)

Galaxy Zoo 1 (Lintott et al. 2011, MNRAS 410, 166), fetched from VizieR table
`J/MNRAS/410/166/galaxies` via the TAP service. "Clean" spiral sample:
`fS=1 AND (pcS>0.8 OR paS>0.8)` — 30,126 galaxies.

Columns: `RAJ2000, DEJ2000` (deg, J2000), `pcS` (clockwise vote fraction),
`paS` (anticlockwise vote fraction), `N` (votes). Spin sign is taken as +1
(clockwise) if `pcS>paS`, else -1.

Caveats: (1) SDSS footprint (northern Galactic cap), NOT all-sky — the dipole
fit is mask-limited. (2) GZ1 handedness carries a known human *perception bias*
(Land et al. 2008): the raw CW/CCW asymmetry is largely a systematic, debiased
by mirrored-image votes. This catalogue is the cautionary de-confounding case,
not a clean cosmological probe.

Expected schema for any drop-in spin catalogue used by `cartasis_sims.galaxy_spins`:
`ra_deg, dec_deg, spin` with `spin in {+1 (CW), -1 (CCW)}` — or GZ1-style
`RAJ2000, DEJ2000, pcS, paS`.

To get an all-sky probe, replace with a full-sky spin catalogue (e.g. Shamir's
DESI Legacy Survey tables) under the same schema.

## gz1_raw_debiased.csv  (REAL DATA)

The same GZ1 clean-spiral cut as `gz1_clean.csv`, but pulling the *debiased*
spin-vote fractions (`pcSm, paSm`) alongside the raw ones (`pcS, paS`). For the
clean sample the debiasing does not flip the spin sign of any galaxy, so the
two derivations give the *same* per-galaxy handedness — GZ1's known asymmetry is
a vote-level perception bias (Land et al. 2008), not removed by these columns.

## iye2019.csv  (REAL DATA — independent method)

Iye, Tadaki & Fukumoto (2019, ApJ 886, 133), "The Spin Parity of Spiral
Galaxies. II", fetched from VizieR tables `J/ApJ/886/133/figure{8,9,10,11}` via
TAP and merged (deduplicated by ID): 530 spirals with expert-determined **S/Z
winding** (chirality) and sky positions. Columns: `ra_deg, dec_deg, spin`
(`+1` = S-wise, `-1` = Z-wise).

This is a genuinely *independent* handedness determination from GZ1 — a small,
expert-curated, broader-sky sample (this is the reanalysis that found the spiral
handedness consistent with isotropy on large scales), versus GZ1's large,
crowd-sourced, northern SDSS sample with its perception bias. Comparing the two
exposes how method and footprint, not the cosmos, drive the inferred spin axis.
Used by `figures/scripts/ch06_spin_datasets.py`.
