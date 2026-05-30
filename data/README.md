# Data

Observational datasets used for comparison with SCT predictions (Chapter 8).
Mirrors the role of `data/` in the Wealth of Agents project: raw extracts plus
provenance, never edited in place.

## Intended subdirectories

```
data/
├── cmb/        # Planck / WMAP power spectra, low-multipole maps ("axis of evil")
├── desi/       # DESI BAO / dark-energy w0-wa constraints
├── jwst/       # JWST galaxy-rotation-direction catalogs (Shamir-style analyses)
└── README.md   # this file
```

## Provenance rule

Each dataset gets a short note recording: source URL, version/date of download,
and any preprocessing applied. Raw files are treated as read-only inputs;
derived products go to `sims/output/` or `figures/data/`.
