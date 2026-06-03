r"""Tier 2: the tensor-to-scalar ratio r and the inflation-vs-bounce fork.

The primordial *tilt* is in primordial.py (Mukhanov--Sasaki through the bounce). This
module adds the second observable, r = P_T/P_S, and frames the clean discriminator
between inflation and the matter/torsion bounce. The point is NOT a single predicted
number -- it is which CONSISTENCY RELATION the data land on:

  * Inflation ties the tensor tilt to r:           n_T = -r/8     (slow-roll).
    Slow-roll plateau models (Starobinsky/Higgs) sit at small r ~ 0.003-0.004 with a
    tiny red n_T; large-field m^2 phi^2 gives r ~ 0.13 (already excluded).

  * The matter bounce ties the tensor tilt to the SCALAR tilt:   n_T = n_s - 1.
    Both spectra come from the same modes leaving the horizon during the contracting
    phase, so their tilts are EQUAL. With n_s = 0.965 this predicts n_T = -0.035 --
    about seventy times more tilt than a plateau inflation model at the same r, and
    decoupled from r entirely.

On r's MAGNITUDE the bounce is honestly not yet sharp: the simplest single-field
matter bounce OVER-produces tensors (r of order unity, excluded by r < 0.036), which
the realistic realisation must cure by enhancing the scalar power (a small scalar
sound speed c_s < 1, or a second field). Once cured, r is a model-dependent output
that generically drops below LiteBIRD's reach. So the near-term fork is:

  - a tensor DETECTION at the plateau value (~0.003) sitting on the inflationary line
    n_T = -r/8 favours inflation;
  - a NULL (r below ~0.001), or a measured n_T ~ n_s-1 far off the inflationary line,
    favours the bounce.

What is robust and computed (primordial.py): n_s ~ 0.965 from a contraction just
softer than matter, and the equal-tilt relation n_T = n_s - 1. What is owed: the
absolute r, which needs the scalar sound-speed / sound-horizon history of the torsion
bounce (Ch.8, Tier 2).
"""

from __future__ import annotations

from . import primordial as pr

# current and forecast tensor sensitivity
BK18_UPPER = 0.036       # r_0.05 < 0.036 (BICEP/Keck 2021, 95% CL)
LITEBIRD_SIGMA_R = 1e-3  # LiteBIRD total target sigma(r) ~ 0.001 (incl. foregrounds)
N_S_OBS = 0.9649         # Planck 2018 scalar tilt
N_S_SIGMA = 0.0042


def starobinsky(n_efold: float = 57.0):
    """(n_s, r) for Starobinsky R^2 / Higgs-like plateau inflation at N e-folds:
    n_s = 1 - 2/N, r = 12/N^2. N~57 gives (0.965, 0.0037)."""
    return 1.0 - 2.0 / n_efold, 12.0 / n_efold**2


def chaotic_phi2(n_efold: float = 57.0):
    """(n_s, r) for m^2 phi^2 large-field inflation: n_s = 1 - 2/N, r = 8/N.
    N~57 gives (0.965, 0.14) -- excluded by r < 0.036."""
    return 1.0 - 2.0 / n_efold, 8.0 / n_efold


def inflation_consistency_nT(r: float) -> float:
    """Inflationary single-field consistency relation: n_T = -r/8."""
    return -r / 8.0


def bounce_tensor_tilt(n_s: float = N_S_OBS) -> float:
    """The matter-bounce equal-tilt relation: n_T = n_s - 1 (scalar and tensor
    spectra share the contracting-phase generation, so their tilts coincide)."""
    return n_s - 1.0


def bounce_tensor_tilt_from_w(w: float) -> float:
    """Tensor tilt of the bounce family from the contraction equation of state w,
    via the Mukhanov--Sasaki evolution (analytic n_T = 3 - |2p-1|, p = 2/(1+3w))."""
    return pr.analytic_tilt(pr.p_from_w(w))


def w_for_observed_tilt(n_s: float = N_S_OBS) -> float:
    """Contraction w reproducing the observed scalar tilt: n_T = n_s-1 = 12 w near
    matter domination, so w = (n_s - 1)/12."""
    return (n_s - 1.0) / 12.0


def discriminator_gap(n_s: float = N_S_OBS, r_inflation: float = None) -> dict:
    """The tensor-tilt gap between the two hypotheses at the observed n_s:
    bounce n_T = n_s-1, inflation n_T = -r/8 (Starobinsky r if r not given)."""
    if r_inflation is None:
        _, r_inflation = starobinsky()
    nT_bounce = bounce_tensor_tilt(n_s)
    nT_infl = inflation_consistency_nT(r_inflation)
    return dict(nT_bounce=nT_bounce, nT_inflation=nT_infl,
                gap=nT_infl - nT_bounce)
