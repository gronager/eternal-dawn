# DESI DR2 dark-energy constraints

CPL (w0waCDM) best-fit values from DESI DR2 BAO (March 2025),
arXiv:2503.14738 (DESI DR2 Results II) and arXiv:2503.14743 (extended analysis).

`desi_dr2_w0wa.csv` columns: dataset, w0, w0_err, wa, wa_err, sigma_LCDM.

CAVEAT (VERIFY): central values/errors below are transcribed from memory of the
DR2 release and flagged for verification against the published tables before any
quotable use. The *qualitative* result is robust and stated by the collaboration:
all combinations favour w0 > -1 and wa < 0 (dynamical DE), 2.8-4.2 sigma over LCDM,
with w(z) phantom (< -1) in the past crossing to > -1 today.

The framework's prediction (Ch. 6): dark energy = response to the parent black
hole's accretion, so w evolves; w(z) ~ -1 + (something) d ln(Mdot_parent)/d ln a.
A still-growing parent gives w0 > -1 with w increasing toward the past being more
negative -> wa < 0. This matches the DESI sign pattern.
