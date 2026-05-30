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
