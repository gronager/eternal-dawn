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


def effective_mass_table(R, T, W):
    """Diagnostic: the effective potential V_eff(R,T) = ln[W(R,T)/W(R,T+1)] as a function of T,
    for each R. This is the raw object whose T-plateau IS the static potential V(R): at small T
    excited states inflate V_eff, then it falls to a plateau, then noise blows it up. Printing
    this table tells you (a) where the plateau is -- the right fit window -- and (b) whether you
    have a plateau at all (no plateau within the signal window => you need smearing).
    Returns {R: (T_mid_array, Veff_array)} with T_mid the lower T of each adjacent pair."""
    R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    out = {}
    for r in np.unique(R):
        m = R == r
        Tr, Wr = T[m], W[m]
        order = np.argsort(Tr)
        Tr, Wr = Tr[order], Wr[order]
        with np.errstate(divide="ignore", invalid="ignore"):
            veff = np.log(Wr[:-1] / Wr[1:])
        out[int(r)] = (Tr[:-1].astype(float), veff)
    return out


def sommer_scale(alpha, sigma, ref=1.65):
    """The Sommer scale r0, defined by r0^2 dV/dr |_{r0} = ref (ref=1.65 conventional). For the
    Cornell form dV/dr = alpha/r^2 + sigma, this gives r0 = sqrt((ref - alpha)/sigma)."""
    if sigma <= 0:
        return float("nan")
    val = (ref - alpha) / sigma
    return float(np.sqrt(val)) if val > 0 else float("nan")


def _potential_and_fit(R, T, W, tmin, tmax):
    """Average W over the (already-selected) rows per (R,T), extract V(R) from the plateau window,
    keep the leading monotonic-rising R range, and Cornell-fit. Returns (fit_dict, R_used) or
    (None, None) if fewer than 3 clean points. Shared by the point estimate and the jackknife."""
    R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    Ru, Tu, Wu = [], [], []
    for (r, t) in sorted(set(zip(R.astype(int), T.astype(int)))):
        m = (R == r) & (T == t)
        Ru.append(r); Tu.append(t); Wu.append(W[m].mean())
    Rr, Vr = effective_potential(np.array(Ru), np.array(Tu), np.array(Wu), t_window=(tmin, tmax))
    good = np.isfinite(Vr) & (Vr > 0)
    keep, last = [], -np.inf
    for i in range(len(Rr)):
        if not good[i]:
            break
        if Vr[i] >= last - 1e-9:
            keep.append(i); last = Vr[i]
        else:
            break
    keep = np.array(keep, dtype=int)
    if keep.size < 3:
        return None, None
    return static_potential_cornell(Rr[keep], Vr[keep]), Rr[keep]


def string_tension_jackknife(cfg, R, T, W, tmin=2, tmax=4):
    """String tension with a delete-1 jackknife error over configurations. `cfg` is the per-row
    configuration id; (R, T, W) are the timelike Wilson loops from every config. The full sample
    gives the central sigma/alpha; deleting one config at a time and refitting gives the jackknife
    error sqrt((n-1)/n * sum (x_i - mean)^2). Error bars are what separate a real measurement
    (sigma small relative to its error) from a noise-dominated fit (error >= the value, e.g. the
    unphysical negative alpha you get from a handful of correlated configs)."""
    cfg = np.asarray(cfg); R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    cfgs = np.unique(cfg)
    n = len(cfgs)
    full_fit, R_used = _potential_and_fit(R, T, W, tmin, tmax)
    if full_fit is None or n < 2:
        return {"sigma": float("nan"), "sigma_err": float("nan"), "alpha": float("nan"),
                "alpha_err": float("nan"), "n_cfg": int(n), "R_used": []}
    js, ja = [], []
    for c in cfgs:
        m = cfg != c
        fit, _ = _potential_and_fit(R[m], T[m], W[m], tmin, tmax)
        if fit is not None:
            js.append(fit["sigma"]); ja.append(fit["alpha"])
    js, ja = np.array(js), np.array(ja)

    def jk_err(vals):
        if vals.size < 2:
            return float("nan")
        return float(np.sqrt((vals.size - 1) / vals.size * np.sum((vals - vals.mean()) ** 2)))

    return {"sigma": float(full_fit["sigma"]), "sigma_err": jk_err(js),
            "alpha": float(full_fit["alpha"]), "alpha_err": jk_err(ja),
            "r0_sommer": float(full_fit["r0_sommer"]), "n_cfg": int(n),
            "R_used": [int(x) for x in R_used]}


def beta_from_plaquette(betas, plaqs, target=0.5937):
    """Map the bare coupling to a target lattice spacing via the plaquette: given equilibrium
    average plaquettes <P>(beta) measured at a handful of input betas, find the beta whose
    equilibrium plaquette equals `target`. Useful for landing on a chosen reference point. Grid's
    `beta` is the standard 6/g^2 (verified: input 6.0 -> <P>=0.5937, 6.2 -> 0.6136, both matching
    the literature), so this is calibration to a spacing, not a convention fix. Reference SU(3)
    Wilson plaquettes: beta=5.7 -> <P>~0.549, 5.8 -> ~0.567, 6.0 -> ~0.5937, 6.2 -> ~0.6136.

    Monotonic in beta, so we interpolate (linearly in beta vs <P>) and, if `target` lies outside
    the scanned range, linearly extrapolate from the two nearest points. Returns the calibrated
    beta and a flag for whether it was interpolated (in-range) or extrapolated."""
    betas = np.asarray(betas, dtype=float)
    plaqs = np.asarray(plaqs, dtype=float)
    order = np.argsort(plaqs)                      # plaquette increases monotonically with beta
    betas, plaqs = betas[order], plaqs[order]
    inside = plaqs.min() <= target <= plaqs.max()
    beta_star = float(np.interp(target, plaqs, betas))
    if not inside:                                 # np.interp clamps; extrapolate from the nearest edge instead
        if target > plaqs.max():
            (p0, p1), (b0, b1) = plaqs[-2:], betas[-2:]
        else:
            (p0, p1), (b0, b1) = plaqs[:2], betas[:2]
        beta_star = float(b0 + (b1 - b0) * (target - p0) / (p1 - p0))
    return {"beta": beta_star, "target_plaq": float(target), "interpolated": bool(inside)}


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
