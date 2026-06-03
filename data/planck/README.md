# Planck 2018 CMB TT power spectrum (binned)

`COM_PowerSpect_CMB-TT-binned_R3.01.txt` — the Planck 2018 (PR3) binned
temperature (TT) angular power spectrum, as released by the Planck Collaboration
through the Planck Legacy Archive.

- **Source:** ESA Planck Legacy Archive, product `COM_PowerSpect_CMB-TT-binned_R3.01.txt`
  (`https://pla.esac.esa.int/`).
- **Columns:** `ell`, `Dl` [μK²], `-dDl`, `+dDl` (lower/upper 1σ), `BestFit`
  (the Planck 2018 base-ΛCDM best-fit `Dl`), where `Dl = ell(ell+1)C_ell/2π`.
- **Coverage:** 83 band powers, ℓ ≈ 48 … 2500.

Used by `figures/scripts/ch09_peak_heights.py` to compare the Eternal Dawn
full-Boltzmann (CAMB) prediction against the *real* measured spectrum and to
read off the acoustic-peak heights. The framework runs the standard post-bounce
recombination, so the comparison is a genuine consistency test, not a re-fit: see
`chapters/08-simulations.tex` §"Tier 2".

This is observational data, redistributed under the Planck Collaboration's terms;
cite Planck 2018 results (Aghanim et al. 2020) when using it.
