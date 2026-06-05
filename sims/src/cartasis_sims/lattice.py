r"""Lattice analysis extractors for the Eternal Dawn strong-sector targets.

The lattice CONFIGURATIONS are generated on the GPU (Grid/QUDA on the GH200; see the
`lattice/` directory). This module is the laptop/CPU-testable ANALYSIS half: it turns the
measurement outputs into the three numbers the programme needs, and each extractor ships with
a synthetic-data test that recovers a known input, so the analysis logic is validated
independently of any GPU run.

The three targets (Appendix B, run order cheapest-first):
  * sigma  -- string tension / confinement, from the static potential V(r)=c-alpha/r+sigma r
              (target L1-L2; pure gauge; the pipeline validator).
  * T_c    -- the deconfinement / chiral-restoration temperature, from the peak of the Polyakov
              susceptibility (the lattice version of 'the condensate melts' -- genesis Stage 3).
  * gamma_m-- the mass anomalous dimension, from the Dirac MODE NUMBER nu(M) ~ M^{4/(1+gamma)}
              (Giusti-Luscher), the GATE measurement that decides whether the candidate walks.

Plus the gradient-flow scale w0 for setting the lattice spacing. References: the static
potential (Cornell/Sommer r0), the Polyakov deconfinement transition, the Giusti-Luscher mode
number, and the Wilson/gradient flow are all standard lattice tools.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 1. String tension from the static potential (target L1-L2)
# ---------------------------------------------------------------------------
def static_potential_cornell(r, V, V_err=None):
    """Fit the static quark potential to the Cornell form V(r) = c - alpha/r + sigma*r and
    return the string tension sigma (and the Coulomb coefficient alpha, constant c, and the
    Sommer scale r0). A nonzero sigma is confinement (the area law); sigma -> 0 is screening.

    Inputs in lattice units; sigma comes out in 1/a^2. Pass V_err for weighted least squares."""
    r = np.asarray(r, dtype=float)
    V = np.asarray(V, dtype=float)
    # linear model in the basis {1, 1/r, r}: V = c*1 + (-alpha)*(1/r) + sigma*r
    A = np.column_stack([np.ones_like(r), 1.0 / r, r])
    if V_err is not None:
        w = 1.0 / np.asarray(V_err, dtype=float)
        A, y = A * w[:, None], V * w
    else:
        y = V
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    c, neg_alpha, sigma = coef
    alpha = -neg_alpha
    return {"sigma": float(sigma), "alpha": float(alpha), "c": float(c),
            "r0_sommer": sommer_scale(alpha, sigma)}


def potential_from_loops_fit(R, T, W, tmin=1, tmax=4):
    """V(R) as the slope of -ln W(R,T) vs T over [tmin, tmax], fit per R. More robust than a
    single effective-mass ratio when statistics are limited (it uses all T in the window and
    is less sensitive to one noisy point). Returns (R, V(R)) for static_potential_cornell.
    Without ground-state smearing the small-T slope over-estimates V (excited-state
    contamination), but it still rises with R -- enough to see confinement and sign(sigma)."""
    R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    Rs = np.unique(R)
    V = []
    for r in Rs:
        m = (R == r) & (T >= tmin) & (T <= tmax) & (W > 0)
        Tr, lw = T[m], -np.log(W[m])
        if Tr.size >= 2 and np.unique(Tr).size >= 2:
            V.append(float(np.polyfit(Tr, lw, 1)[0]))
        else:
            V.append(float("nan"))
    return Rs.astype(float), np.array(V)


def effective_potential(R, T, W, t_window=None):
    """The static potential V(R) from timelike Wilson loops W(R,T), via the effective-mass
    ratio V(R) = ln[W(R,T)/W(R,T+1)], averaged over a plateau window in T (default: the upper
    half of the available T, where excited-state contamination has decayed). Returns (R, V(R))
    ready for static_potential_cornell."""
    R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    Rs = np.unique(R)
    Vr = []
    for r in Rs:
        m = R == r
        Tr, Wr = T[m], W[m]
        order = np.argsort(Tr)
        Tr, Wr = Tr[order], Wr[order]
        vt = np.log(Wr[:-1] / Wr[1:])              # effective V at each T
        Tmid = Tr[:-1]
        if t_window is None:
            sel = Tmid >= max(1, int(Tr.max()) // 2)
        else:
            sel = (Tmid >= t_window[0]) & (Tmid <= t_window[1])
        vt_sel = vt[sel & np.isfinite(vt)]
        Vr.append(float(np.mean(vt_sel)) if vt_sel.size else float("nan"))
    return Rs.astype(float), np.array(Vr)


def sommer_scale(alpha, sigma, ref=1.65):
    """The Sommer scale r0, defined by r0^2 dV/dr |_{r0} = ref (ref=1.65 conventional). For the
    Cornell form dV/dr = alpha/r^2 + sigma, this gives r0 = sqrt((ref - alpha)/sigma)."""
    if sigma <= 0:
        return float("nan")
    val = (ref - alpha) / sigma
    return float(np.sqrt(val)) if val > 0 else float("nan")


# ---------------------------------------------------------------------------
# 2. Deconfinement temperature from the Polyakov susceptibility (the melting)
# ---------------------------------------------------------------------------
def deconfinement_beta_c(beta, chi_L):
    """The critical coupling beta_c at the deconfinement transition: the peak of the Polyakov-
    loop susceptibility chi_L(beta). Located by a parabola fit around the discrete maximum, so
    the result is not pinned to the scan grid. Returns beta_c and the peak height."""
    beta = np.asarray(beta, dtype=float)
    chi = np.asarray(chi_L, dtype=float)
    i = int(np.argmax(chi))
    if 0 < i < len(beta) - 1:
        x = beta[i - 1:i + 2]
        y = chi[i - 1:i + 2]
        a, b, c = np.polyfit(x, y, 2)
        beta_c = -b / (2 * a)
        peak = a * beta_c**2 + b * beta_c + c
    else:
        beta_c, peak = beta[i], chi[i]
    return {"beta_c": float(beta_c), "chi_peak": float(peak)}


# ---------------------------------------------------------------------------
# 3. Mass anomalous dimension from the Dirac mode number (the GATE)
# ---------------------------------------------------------------------------
def anomalous_dimension_from_mode_number(M, nu, window=None):
    """The mass anomalous dimension gamma_m from the Dirac mode number nu(M) (the number of
    eigenvalues of the massless Dirac operator below M). In the near-conformal scaling regime
    nu(M) ~ M^{4/(1+gamma_m)} (Giusti-Luscher), so the local log-log slope
        s = d ln nu / d ln M = 4/(1+gamma_m)  ->  gamma_m = 4/s - 1.
    `window` = (M_lo, M_hi) restricts the fit to the scaling window (recommended). Returns
    gamma_m and the fitted slope."""
    M = np.asarray(M, dtype=float)
    nu = np.asarray(nu, dtype=float)
    good = (M > 0) & (nu > 0)
    if window is not None:
        good &= (M >= window[0]) & (M <= window[1])
    lnM, lnnu = np.log(M[good]), np.log(nu[good])
    slope, _ = np.polyfit(lnM, lnnu, 1)
    gamma_m = 4.0 / slope - 1.0
    return {"gamma_m": float(gamma_m), "slope": float(slope)}


# ---------------------------------------------------------------------------
# 4. Gradient-flow scale w0 (sets the lattice spacing)
# ---------------------------------------------------------------------------
def gradient_flow_w0(t, t2E, ref=0.3):
    """The gradient-flow scale w0, defined by W(t) = t d/dt [t^2 <E>(t)] = ref at t = w0^2
    (ref=0.3 conventional). Returns w0 in lattice units (so a = w0_phys / w0_lat)."""
    t = np.asarray(t, dtype=float)
    f = np.asarray(t2E, dtype=float)
    W = t * np.gradient(f, t)
    # first crossing of W = ref, linearly interpolated
    above = W >= ref
    idx = np.where(above[1:] & ~above[:-1])[0]
    if len(idx) == 0:
        idx = np.where(above[:-1] & ~above[1:])[0]   # falling crossing fallback
    if len(idx) == 0:
        return float("nan")
    i = int(idx[0])
    t_star = t[i] + (ref - W[i]) * (t[i + 1] - t[i]) / (W[i + 1] - W[i])
    return float(np.sqrt(t_star))


CANDIDATE_THEORY = {
    "gauge_group": "SU(3)",
    "fermions": "N_f = 2 Dirac in the sextet (2-index symmetric) representation",
    "why": "leading near-conformal/walking candidate with a light flavour-singlet scalar "
           "(the composite-Higgs/dilaton our Part III needs); SU(3) matches the colour the "
           "Pauli label forces; substantial existing lattice literature to validate against.",
    "gate_observable": "gamma_m via the Dirac mode number",
    "expected_gamma_m": "order 0.3 if walking (literature, debated); the sign that decides it.",
}
