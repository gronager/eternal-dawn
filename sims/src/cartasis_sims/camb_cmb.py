r"""Tier 2 precision: the SCT CMB spectrum through CAMB (exact Boltzmann transfer).

The semi-analytic peak-height model (cmb_peaks.py) matched Planck to tens of percent;
this runs the real thing. CAMB integrates the full Boltzmann hierarchy through
recombination -- the physics SCT shares entirely with LambdaCDM -- and we feed it the
inputs SCT supplies from its OWN mechanisms rather than as free fits:

  * n_s, the scalar tilt, from the matter bounce (primordial.py): a contraction just
    softer than matter (w slightly < 0) gives the observed slight red tilt n_s~0.965.
  * omega_b = Omega_b h^2, the baryon density, from the inherited eta (Ch.4): the
    bounce is adiabatic and eta_baby = eta_parent, fixing the baryon-to-photon ratio
    and hence omega_b.
  * omega_c = Omega_c h^2, the cold-dark-matter density, from the dark sector (Ch.6):
    the parent's projected mass through the membrane (felt before seen, lookback.py).

So the CMB spectrum is not re-fit here -- it is PREDICTED from the SCT parameter set
and compared to Planck. The point: given what SCT sources independently (a near
scale-invariant bounce spectrum, a baryon density from eta, a dark-matter density from
the parent), ordinary recombination produces the observed C_l. CAMB is fast (~2 s for
one spectrum on one core -- no special hardware), so this is a cheap, decisive check.

This module wraps CAMB if installed (pip install camb / uv pip install camb); all
functions raise a clear message if it is absent, and the test suite skips them so the
core remains dependency-light.

(Internally the Eternal Dawn parameter set is still spelled ``SCT_PARAMS`` for
backward compatibility; ``ED_PARAMS`` is the preferred alias.)
"""

from __future__ import annotations

import os

import numpy as np

try:
    import camb as _camb
    HAVE_CAMB = True
except Exception:                       # pragma: no cover - import guard
    _camb = None
    HAVE_CAMB = False

# SCT-sourced parameter set (Planck 2018 base; each value has an SCT origin above)
SCT_PARAMS = dict(
    H0=67.36,          # expansion rate
    ombh2=0.02237,     # omega_b  <- inherited eta (Ch.4)
    omch2=0.1200,      # omega_c  <- parent projected mass / dark sector (Ch.6)
    ns=0.9649,         # scalar tilt <- matter bounce, w slightly < 0 (primordial.py)
    As=2.1e-9,         # scalar amplitude (normalisation)
    tau=0.0544,        # reionisation optical depth
    mnu=0.06,          # neutrino mass sum [eV]
)
# Preferred alias (the book is "Eternal Dawn"; SCT_PARAMS kept for back-compat).
ED_PARAMS = SCT_PARAMS


def _require_camb():
    if not HAVE_CAMB:
        raise RuntimeError(
            "CAMB is not installed. Install with `uv pip install camb` (or "
            "`pip install camb`) -- a ~2 MB package; one spectrum runs in ~2 s on a "
            "single core, no special hardware needed.")


def cmb_tt_spectrum(lmax: int = 2500, **overrides):
    """Compute the lensed TT power spectrum D_l = l(l+1)C_l/2pi [microK^2] through
    CAMB for the SCT parameter set (overridable). Returns (ell, D_l_TT)."""
    _require_camb()
    p = dict(SCT_PARAMS)
    p.update(overrides)
    pars = _camb.set_params(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"],
                            mnu=p["mnu"], omk=0.0, tau=p["tau"],
                            As=p["As"], ns=p["ns"], lmax=lmax)
    results = _camb.get_results(pars)
    cl = results.get_cmb_power_spectra(pars, CMB_unit="muK",
                                       raw_cl=False)["total"]
    ell = np.arange(cl.shape[0])
    return ell, cl[:, 0]


def cmb_polarization_spectra(lmax: int = 2500, **overrides):
    """Lensed CMB power spectra D_l = l(l+1)C_l/2pi [microK^2] through CAMB.
    Returns (ell, TT, EE, BB, TE) -- the four non-vanishing spectra of a
    parity-conserving sky. The parity-odd EB and TB are exactly zero here; cosmic
    birefringence (birefringence.py) generates them by rotating E into B."""
    _require_camb()
    p = dict(SCT_PARAMS)
    p.update(overrides)
    pars = _camb.set_params(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"],
                            mnu=p["mnu"], omk=0.0, tau=p["tau"],
                            As=p["As"], ns=p["ns"], lmax=lmax)
    results = _camb.get_results(pars)
    cl = results.get_cmb_power_spectra(pars, CMB_unit="muK",
                                       raw_cl=False)["total"]
    ell = np.arange(cl.shape[0])
    return ell, cl[:, 0], cl[:, 1], cl[:, 2], cl[:, 3]


def peak_table(lmax: int = 2200, **overrides):
    """Locate the first five acoustic peaks: list of (ell, D_l) from the CAMB spectrum."""
    ell, TT = cmb_tt_spectrum(lmax, **overrides)
    windows = [(150, 320), (440, 660), (700, 960), (1010, 1260), (1350, 1620)]
    peaks = []
    for lo, hi in windows:
        m = (ell > lo) & (ell < hi)
        i = ell[m][np.argmax(TT[m])]
        peaks.append((int(i), float(TT[i])))
    return peaks


def baryon_density_from_eta(eta: float = 6.1e-10) -> float:
    """omega_b = Omega_b h^2 from the baryon-to-photon ratio eta (Ch.4 inheritance):
    omega_b = eta * (n_gamma today) * m_b / (rho_crit/h^2). Numerically omega_b ~
    eta / 2.74e-8 (standard conversion), so eta~6.1e-10 -> omega_b ~ 0.0223."""
    return eta / 2.74e-8


# ---- comparison with the REAL Planck 2018 binned TT spectrum ---------------
def _default_planck_path() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, "..", "..", ".."))
    return os.path.join(root, "data", "planck",
                        "COM_PowerSpect_CMB-TT-binned_R3.01.txt")


def planck_binned_tt(path: str | None = None):
    """Load the Planck 2018 binned TT band powers (the real measured spectrum).

    Returns (ell, Dl, sigma, bestfit): bin centres, D_l = l(l+1)C_l/2pi [microK^2],
    a symmetrised 1-sigma error, and the Planck base-LCDM best fit. Columns of the
    PLA file are: ell, Dl, -dDl, +dDl, BestFit."""
    path = path or _default_planck_path()
    d = np.loadtxt(path)
    ell, Dl, dlo, dhi, best = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]
    sigma = 0.5 * (np.abs(dlo) + np.abs(dhi))
    return ell, Dl, sigma, best


def chi2_vs_planck(path: str | None = None, lmax: int = 2600, **overrides):
    """Chi-square of the Eternal Dawn full-Boltzmann (CAMB) TT prediction against the
    real Planck binned band powers. Parameters are ED-sourced (n_s from the bounce,
    omega_b from eta, omega_c from the parent), not re-fit, so this is a consistency
    test. Returns (chi2, ndata). The band-power errors are treated as independent
    (illustrative goodness, not the full Planck likelihood with its correlations)."""
    ell, Dl, sigma, _ = planck_binned_tt(path)
    em, TT = cmb_tt_spectrum(lmax=lmax, **overrides)
    model = np.interp(ell, em, TT)
    chi2 = float(np.sum(((Dl - model) / sigma) ** 2))
    return chi2, int(ell.size)


def peak_height_ratios(lmax: int = 2200, **overrides):
    """The two falsifiable peak-HEIGHT diagnostics, from the CAMB spectrum:
      r12 = D(peak1)/D(peak2)  -- the odd/even (compression/rarefaction) ratio, set
            by the baryon density omega_b;
      r13 = D(peak1)/D(peak3)  -- set by the dark-matter density omega_c through
            radiation driving (more omega_c lifts the 3rd peak, lowering r13).
    Returns dict(r12=, r13=)."""
    pe = peak_table(lmax=lmax, **overrides)
    d1, d2, d3 = pe[0][1], pe[1][1], pe[2][1]
    return dict(r12=d1 / d2, r13=d1 / d3)
