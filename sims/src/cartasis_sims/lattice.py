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


def _potential_and_fit(R, T, W, tmin, tmax, rmin=1, rmax=None):
    """Average W over the (already-selected) rows per (R,T), extract V(R) from the plateau window,
    restrict to the fit window [rmin, rmax], keep the leading monotonic-rising R range, and
    Cornell-fit. Returns (fit_dict, R_used) or (None, None) if fewer than 3 clean points. Shared
    by the point estimate and the jackknife. Capping rmax matters: the largest loops (R ~ L/2)
    are the noisiest, and letting them flip in and out of the kept range across jackknife samples
    inflates the error -- excluding them is what stabilises the fit."""
    R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    Ru, Tu, Wu = [], [], []
    for (r, t) in sorted(set(zip(R.astype(int), T.astype(int)))):
        m = (R == r) & (T == t)
        Ru.append(r); Tu.append(t); Wu.append(W[m].mean())
    Rr, Vr = effective_potential(np.array(Ru), np.array(Tu), np.array(Wu), t_window=(tmin, tmax))
    if rmax is None:
        rmax = np.inf
    good = np.isfinite(Vr) & (Vr > 0) & (Rr >= rmin) & (Rr <= rmax)
    keep, last = [], -np.inf
    for i in range(len(Rr)):
        if Rr[i] < rmin:                  # skip below the window without breaking the run
            continue
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


def string_tension_jackknife(cfg, R, T, W, tmin=2, tmax=4, rmin=1, rmax=None):
    """String tension with a delete-1 jackknife error over configurations. `cfg` is the per-row
    configuration id; (R, T, W) are the timelike Wilson loops from every config. The full sample
    gives the central sigma/alpha; deleting one config at a time and refitting gives the jackknife
    error sqrt((n-1)/n * sum (x_i - mean)^2). Error bars are what separate a real measurement
    (sigma small relative to its error) from a noise-dominated fit (error >= the value, e.g. the
    unphysical negative alpha you get from a handful of correlated configs)."""
    cfg = np.asarray(cfg); R = np.asarray(R); T = np.asarray(T); W = np.asarray(W, dtype=float)
    cfgs = np.unique(cfg)
    n = len(cfgs)
    full_fit, R_used = _potential_and_fit(R, T, W, tmin, tmax, rmin, rmax)
    if full_fit is None or n < 2:
        return {"sigma": float("nan"), "sigma_err": float("nan"), "alpha": float("nan"),
                "alpha_err": float("nan"), "n_cfg": int(n), "R_used": []}
    js, ja = [], []
    for c in cfgs:
        m = cfg != c
        fit, _ = _potential_and_fit(R[m], T[m], W[m], tmin, tmax, rmin, rmax)
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
def mode_number_from_eigenvalues(eigenvalues, M_grid):
    """nu(M) = number of eigenvalues of the Hermitian Dirac operator with |lambda| < M, as a
    function of the threshold M. The mode number is just the cumulative count of the spectrum;
    this is the bridge from a measured (or analytic) Dirac spectrum to the gamma_m fit. Per
    config -- average nu(M) over configs before fitting. eigenvalues are |lambda| (>=0)."""
    eig = np.sort(np.abs(np.asarray(eigenvalues, dtype=float)))
    M = np.asarray(M_grid, dtype=float)
    return np.searchsorted(eig, M, side="right").astype(float)


def free_wilson_mode_number(L, T, M_grid, m=0.0, r=1.0):
    """Analytic free Wilson-Dirac mode number nu(M) on an L^3 x T lattice -- the exact, gauge-free
    reference that the measurement chain must reproduce. The free spectrum is diagonal in momentum
    space: for lattice momentum p_mu = 2*pi*n_mu/N_mu the Hermitian operator has
        |lambda(p)| = sqrt( sum_mu sin^2 p_mu + W(p)^2 ),   W(p) = m + r * sum_mu (1 - cos p_mu),
    each with multiplicity 12 (4 spin x 3 colour). A free (non-interacting) theory has gamma_m = 0,
    so nu(M) ~ M^4 at small M; feeding this to anomalous_dimension_from_mode_number must return
    gamma_m ~ 0. This is the first rung: validate the nu(M) -> gamma_m pipeline against an exact
    answer before trusting any interacting (Grid) measurement, and long before the sextet."""
    lam = np.sort(_free_wilson_abs_lambda(L, T, m=m, r=r))
    M = np.asarray(M_grid, dtype=float)
    return 12.0 * np.searchsorted(lam, M, side="right").astype(float)   # x12 spin*colour


def _free_wilson_abs_lambda(L, T, m=0.0, r=1.0):
    """The |lambda| of the free Wilson-Dirac operator on an L^3 x T lattice, one entry per lattice
    momentum (spin x colour multiplicity of 12 applied by the caller): the diagonal, gauge-free
    spectrum |lambda(p)| = sqrt(sum_mu sin^2 p_mu + W(p)^2), W(p) = m + r sum_mu (1 - cos p_mu)."""
    axes = [2.0 * np.pi * np.arange(n) / n for n in (L, L, L, T)]
    P = np.meshgrid(*axes, indexing="ij")
    s2 = sum(np.sin(p) ** 2 for p in P)
    W = m + r * sum(1.0 - np.cos(p) for p in P)
    return np.sqrt(s2 + W**2).ravel()


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
# 3b. Mode number the scalable way: stochastic Chebyshev moments (KPM)
# ---------------------------------------------------------------------------
# The Lanczos route (compute the low eigenvalues, count them) needs a high-order Chebyshev filter
# and convergence babysitting. The Kernel Polynomial Method sidesteps both: estimate the Chebyshev
# moments mu_k = Tr T_k(Xtilde) of the rescaled operator Xtilde = (D^dag D)/(xmax/2) - 1 in [-1,1]
# by a stochastic trace (random sources, D^dag D matvecs only -- no eigenvalues, no solves), then
# combine them with the Jackson-damped Chebyshev expansion of the step theta(M^2 - x) to get the
# mode number nu(M). One set of moments yields nu(M) for the whole M grid. This is the division of
# labour with the Grid code: Grid supplies the moments; the functions below do the rest.
def jackson_kernel(N):
    """Jackson damping factors g_0..g_N that suppress the Gibbs oscillations of an order-N
    Chebyshev expansion of a discontinuous function (here the spectral step). Standard KPM kernel
    (Weisse et al., Rev. Mod. Phys. 78, 275)."""
    k = np.arange(N + 1)
    t = np.pi / (N + 1)
    return ((N - k + 1) * np.cos(k * t) + np.sin(k * t) / np.tan(t)) / (N + 1)


def chebyshev_moments_from_eigenvalues(eigenvalues, N, xmax):
    """Exact Chebyshev moments mu_k = sum_i T_k(xtilde_i), xtilde = (lambda^2)/(xmax/2) - 1, for a
    known Dirac spectrum (eigenvalues = |lambda|). This is the noise-free limit of the stochastic
    trace the Grid measurement estimates -- used to validate the moments -> nu(M) pipeline against
    an exact answer before any gauge field."""
    x = np.asarray(eigenvalues, dtype=float) ** 2
    half = xmax / 2.0
    xt = np.clip(x / half - 1.0, -1.0, 1.0)
    mu = np.empty(N + 1)
    Tkm1 = np.ones_like(xt)
    mu[0] = Tkm1.sum()
    if N >= 1:
        Tk = xt.copy()
        mu[1] = Tk.sum()
        for k in range(2, N + 1):
            Tkm1, Tk = Tk, 2.0 * xt * Tk - Tkm1
            mu[k] = Tk.sum()
    return mu


def mode_number_from_chebyshev_moments(moments, M_grid, xmax):
    """nu(M) from Chebyshev moments mu_k of the rescaled operator Xtilde = (D^dag D)/(xmax/2) - 1.
    Combines them with the Jackson-damped Chebyshev expansion of the step theta(M^2 - x):
        nu(M) = sum_k g_k c_k(M) mu_k,
        c_0 = 1 - arccos(Mt)/pi,  c_k = -2 sin(k arccos Mt)/(k pi),  Mt = M^2/(xmax/2) - 1.
    `xmax` must exceed the largest eigenvalue of D^dag D (the spectral bound). Feed the result to
    anomalous_dimension_from_mode_number. This is the cheap post-processing of the Grid moments;
    one moment set covers the entire M grid."""
    mu = np.asarray(moments, dtype=float)
    N = len(mu) - 1
    g = jackson_kernel(N)
    half = xmax / 2.0
    M = np.asarray(M_grid, dtype=float)
    Mt = np.clip(M ** 2 / half - 1.0, -1.0, 1.0)
    ac = np.arccos(Mt)
    k = np.arange(1, N + 1)
    nu = np.empty_like(M)
    for j in range(len(M)):
        c = np.empty(N + 1)
        c[0] = 1.0 - ac[j] / np.pi
        c[1:] = -2.0 * np.sin(k * ac[j]) / (k * np.pi)
        nu[j] = float((g * c * mu).sum())
    return nu


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


def _effective_mass(C):
    """m_eff(t) = ln[C(t)/C(t+1)] for a correlator C(t); length len(C)-1 (index t)."""
    C = np.asarray(C, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(C[:-1] / C[1:])


def baryon_spectrum(raw, T, tmin=None, tmax=None):
    """The pion and nucleon (torsiton) masses from quenched valence correlators (L4 pilot).

    `raw` is the (n_rows, 4) array of rows [cfg, t, C_pi(t), C_N(t)] (measure_baryon output, tagged
    per config). Returns {'pion': ..., 'nucleon': ...}, each a dict with the config-averaged
    effective mass m_eff(t)=ln[C(t)/C(t+1)], and a plateau-window mass with a jackknife-over-configs
    error. Default window [T//8, T//3] avoids small-t excited states and large-t noise; read the
    printed m_eff(t) and narrow it if the plateau sits elsewhere. A clean plateau at m_N > 0 is the
    torsiton's rest mass -- the bound ground state, found non-perturbatively."""
    raw = np.asarray(raw, dtype=float)
    cfg = raw[:, 0].astype(int)
    tt = raw[:, 1].astype(int)
    ts = np.arange(tt.max() + 1)
    configs = np.unique(cfg)
    n = len(configs)
    lo = tmin if tmin is not None else max(1, T // 8)
    hi = tmax if tmax is not None else max(lo + 1, T // 3)
    out = {}
    for name, col in (("pion", 2), ("nucleon", 3)):
        C = raw[:, col]
        Cbar = np.array([C[tt == s].mean() for s in ts])
        meff = _effective_mass(Cbar)
        idx = np.arange(lo, min(hi, len(meff)))
        idx = idx[np.isfinite(meff[idx])]
        res = {"t": ts[:-1], "meff": meff, "tmin": int(lo), "tmax": int(hi), "n_cfg": n}
        if len(idx) < 2:
            res.update(mass=float("nan"), mass_err=float("nan"))
        else:
            mass = float(np.mean(meff[idx]))
            js = []
            for cc in configs:                       # jackknife: drop one config at a time
                keep = cfg != cc
                Cj = np.array([C[keep & (tt == s)].mean() for s in ts])
                mj = _effective_mass(Cj)
                js.append(np.mean(mj[idx]))
            js = np.array(js)
            err = float(np.sqrt((n - 1) / n * np.sum((js - js.mean()) ** 2))) if n > 1 else 0.0
            res.update(mass=mass, mass_err=err)
        out[name] = res
    return out


def _bag_rho(prof, plateau):
    """Config- and plateau-averaged radial scalar density rho(r) from measure_bag_profile rows.

    `prof` is an (n,5) array of [cfg, r2, t, rho_sum, cnt] (shell-summed scalar density and the
    site count of that r^2 shell on that time slice). Returns (r, rho), r=sqrt(r2) in lattice units,
    sorted, with rho the shell-AVERAGED density (sum/count) summed over configs and plateau slices,
    normalised so rho(r_min)=1."""
    lo, hi = plateau
    m = (prof[:, 2] >= lo) & (prof[:, 2] <= hi)
    p = prof[m]
    r2 = np.unique(p[:, 1]).astype(int)
    r, rho = [], []
    for b in r2:
        sel = p[:, 1].astype(int) == b
        cnt = p[sel, 4].sum()
        if cnt <= 0:
            continue
        r.append(np.sqrt(b))
        rho.append(p[sel, 3].sum() / cnt)              # shell average over configs+plateau
    r = np.asarray(r, dtype=float)
    rho = np.asarray(rho, dtype=float)
    order = np.argsort(r)
    r, rho = r[order], rho[order]
    if rho[0] != 0:
        rho = rho / rho[0]                              # normalise to the peak (r->0)
    return r, rho


def _half_density_radius(r, rho):
    """The radius where the (normalised, peak=1) profile falls to 1/2, by linear interpolation in
    the first downward crossing -- the bag's size R0 in lattice units."""
    half = 0.5
    for i in range(1, len(rho)):
        if rho[i] <= half <= rho[i - 1]:
            f = (rho[i - 1] - half) / (rho[i - 1] - rho[i] + 1e-300)
            return float(r[i - 1] + f * (r[i] - r[i - 1]))
    return float(r[-1])                                 # never crosses: bag larger than the box


def bag_profile(prof, T, plateau=(4, 10), r0_over_a=3.166):
    """The torsiton bag shape and the lepton lever (Eternal Dawn, Ch. 11 Generations).

    From the gauge-invariant scalar density rho(r)=Tr[S(x;0)^dag S(x;0)] (measure_bag_profile), the
    config-averaged plateau profile is the dressed-constituent BAG. Its half-density radius R0 (the
    bag size) in units of the Sommer scale r0 gives a dimensionless sharpness

        s_T  ~  R0 / r0 ,

    the size of the mass-giving core relative to the confinement length -- the lever input of
    cartasis_sims.fermion_masses. The lever (harmonic overlap) reproduces the full span only in a
    PRODUCTIVE window s_T in [0.43, 0.70] r0 (it peaks near 0.5 and turns over for sharper cores -- a
    Gaussian-overlap model artifact, so a too-small s_T does NOT mean a bigger hierarchy). The verdict
    flags whether the measured bag lands in that window, and never quotes an out-of-window span as a
    result.

    HONEST on the dictionary: rho(r) is the one-body constituent profile, and the map s_T ~ R0/r0 is a
    first-order identification between the lattice bag and the lever's well-relative core (the absolute
    normalisation would be fixed by matching the Woods-Saxon bag to the overlap model; the genuine
    three-body condensate <N|qbar q(r)|N> is the sequential-source 3-point refinement). What is ROBUST
    here is the MEASURED profile rho(r), R0 in r0 units, and the trend: sharp bag -> large span.

    Returns rho(r), R0 (lattice units), the Fermi-fit wall thickness `a` if it converges, s_T, the
    lever span, jackknife errors over configs, and a verdict."""
    prof = np.asarray(prof, dtype=float)
    from . import fermion_masses as fm

    def lever_span(s):
        # the harmonic overlap lever is NOT monotonic: span rises as the core sharpens only down to a
        # peak near s~0.5, then turns over for s<~0.43 (a model artifact of the Gaussian-core overlap).
        # Trust it in the productive window s in [0.43, 0.70] r0; flag, don't quote-as-fact, outside it.
        sc = float(np.clip(s, 0.30, 1.5))
        lad = fm._ascending_ladder(sc, well=fm.WELL)
        return float(lad[-1] / lad[0])

    r, rho = _bag_rho(prof, plateau)
    R0 = _half_density_radius(r, rho)
    s_T = R0 / r0_over_a
    span = lever_span(s_T)

    # optional Fermi (Woods-Saxon) wall thickness, a secondary shape number
    wall = float("nan")
    try:
        from scipy.optimize import curve_fit
        fermi = lambda rr, R, a: 1.0 / (1.0 + np.exp((rr - R) / a))
        popt, _ = curve_fit(fermi, r, rho, p0=[max(R0, 0.5), 0.5],
                            bounds=([0.1, 0.05], [r.max(), r.max()]), maxfev=20000)
        wall = float(popt[1])
    except Exception:
        pass

    # jackknife over configs for s_T and span
    cfgs = np.unique(prof[:, 0]).astype(int)
    n = len(cfgs)
    js_sT, js_span = [], []
    for c in cfgs:
        rj, rhoj = _bag_rho(prof[prof[:, 0].astype(int) != c], plateau)
        R0j = _half_density_radius(rj, rhoj)
        s = R0j / r0_over_a
        js_sT.append(s)
        js_span.append(lever_span(s))
    js_sT, js_span = np.array(js_sT), np.array(js_span)
    if n > 1:
        sT_err = float(np.sqrt((n - 1) / n * np.sum((js_sT - js_sT.mean()) ** 2)))
        span_err = float(np.sqrt((n - 1) / n * np.sum((js_span - js_span.mean()) ** 2)))
    else:
        sT_err = span_err = float("nan")

    lo_p, hi_p = 0.43, 0.70                            # the lever's trustworthy (productive) window
    under_resolved = R0 < 1.5                          # half-density inside ~1 spacing: barely resolved
    note = (" NOTE: bag barely resolved (R0 < 1.5a to half-density) -- heavy-mass / UV-limited;"
            " read the chiral/continuum trend, not this single point.") if under_resolved else ""

    if s_T < lo_p:
        verdict = (f"bag SMALL/SHARP (R0={R0:.2f}a = {s_T:.2f} r0), BELOW the lever's productive window "
                   f"[{lo_p},{hi_p}] r0. The harmonic lever turns over here, so its span read-out "
                   f"({span:.0f}) is a MODEL ARTIFACT, not a result -- do not read it as the hierarchy. "
                   f"At this heavy pion the constituent is tightly localized; lighter quarks grow the "
                   f"bag toward the productive window. Chase the mass dependence and the 3-pt "
                   f"condensate.{note}")
    elif s_T <= hi_p:
        ratio = span / 3477.0
        tag = ("consistent with DERIVED" if 1.0 / 3 < ratio < 3
               else ("over-shoots" if ratio >= 3 else "under-shoots"))
        verdict = (f"bag in the PRODUCTIVE window (s_T={s_T:.2f} r0): lever span {span:.0f} vs observed "
                   f"3477 -- {tag} (factor {max(ratio,1/ratio):.1f}).{note}")
    else:
        verdict = (f"bag BROAD (s_T={s_T:.2f} r0): span {span:.0f} << 3477 -- the minimal bag is too "
                   f"soft; the mechanism needs sharper IR dynamics (or it fails).{note}")

    return {
        "r": r, "rho": rho, "R0": R0, "wall": wall, "s_T": s_T, "span": span,
        "s_T_err": sT_err, "span_err": span_err, "n_cfg": n, "verdict": verdict,
        "r0_over_a": r0_over_a, "plateau": plateau,
        "under_resolved": bool(under_resolved), "productive_window": (lo_p, hi_p),
    }


def correlator_mass(raw, T, tmin=None, tmax=None):
    """Plateau effective mass of a single correlator C(t), with a jackknife-over-configs error.

    `raw` is (n,3) rows [cfg, t, C(t)] (e.g. the pion C_pi=sum_x rho from measure_bag_profile). Returns
    {mass, mass_err, meff, t, n_cfg} from m_eff(t)=ln[C(t)/C(t+1)] averaged over the window [tmin,tmax]
    (defaults [T//8, T//3]). The lean, one-correlator sibling of baryon_spectrum, reused to read m_pi
    at each valence mass in the bag-profile chiral scan."""
    raw = np.asarray(raw, dtype=float)
    cfg = raw[:, 0].astype(int)
    tt = raw[:, 1].astype(int)
    C = raw[:, 2]
    ts = np.arange(tt.max() + 1)
    configs = np.unique(cfg)
    n = len(configs)
    lo = tmin if tmin is not None else max(1, T // 8)
    hi = tmax if tmax is not None else max(lo + 1, T // 3)
    Cbar = np.array([C[tt == s].mean() for s in ts])
    meff = _effective_mass(Cbar)
    idx = np.arange(lo, min(hi, len(meff)))
    idx = idx[np.isfinite(meff[idx])]
    if len(idx) < 1:
        return {"mass": float("nan"), "mass_err": float("nan"), "meff": meff, "t": ts[:-1], "n_cfg": n}
    mass = float(np.mean(meff[idx]))
    js = []
    for cc in configs:
        keep = cfg != cc
        Cj = np.array([C[keep & (tt == s)].mean() for s in ts])
        js.append(np.mean(_effective_mass(Cj)[idx]))
    js = np.array(js)
    err = float(np.sqrt((n - 1) / n * np.sum((js - js.mean()) ** 2))) if n > 1 else 0.0
    return {"mass": mass, "mass_err": err, "meff": meff, "t": ts[:-1], "n_cfg": n}


# ---------------------------------------------------------------------------
# 5. The pion decay constant f_pi and the substrate scale v (target L2)
# ---------------------------------------------------------------------------
def _cosh_fit_pp(Cbar, m_q, T, tmin=None, tmax=None):
    """Two-sided cosh fit of a zero-momentum pseudoscalar-pseudoscalar correlator C_PP(t) to
        C_PP(t) = (G_PP / (2 m_pi)) * (e^{-m_pi t} + e^{-m_pi (T-t)}),
    returning (m_pi*a, G_PP, f_pi*a). The mass comes from the effective-mass plateau of the cosh
    (acosh-based, periodic) in the fit window; the amplitude G_PP is then the least-squares overlap
    of the (fixed-shape) cosh ansatz with the data over the same window. f_pi*a follows from the
    PCAC / pseudoscalar-density relation (see decay_constant)."""
    C = np.asarray(Cbar, dtype=float)
    if tmin is None:
        tmin = max(1, T // 8)
    if tmax is None:
        tmax = max(tmin + 1, T // 2 - 1)
    tt = np.arange(len(C))
    # periodic effective mass from the cosh ratio: cosh(m(t-T/2))/cosh(m(t+1-T/2)) = C(t)/C(t+1).
    # m_eff(t) solves R = cosh(m(t-T/2))/cosh(m(t+1-T/2)); use the symmetric-cosh acosh estimator
    #   m_eff(t) = arccosh( (C(t-1)+C(t+1)) / (2 C(t)) ),
    # which is exact for a single periodic cosh and avoids choosing the T/2 origin explicitly.
    meff = np.full(len(C), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        for t in range(1, len(C) - 1):
            num = C[t - 1] + C[t + 1]
            den = 2.0 * C[t]
            if den != 0 and num / den >= 1.0:
                meff[t] = np.arccosh(num / den)
    idx = np.arange(tmin, min(tmax, len(C) - 1))
    idx = idx[np.isfinite(meff[idx])]
    if len(idx) < 1:
        return float("nan"), float("nan"), float("nan")
    m_pi = float(np.mean(meff[idx]))
    # amplitude: fit C_PP(t) = (G/(2 m)) * shape(t), shape = e^{-m t}+e^{-m(T-t)}, by least squares
    shape = np.exp(-m_pi * tt) + np.exp(-m_pi * (T - tt))
    sel = idx
    A = (shape[sel] / (2.0 * m_pi))
    denom = float(np.dot(A, A))
    G_PP = float(np.dot(A, C[sel]) / denom) if denom > 0 else float("nan")
    f_pi = _f_pi_from_pcac(m_q, m_pi, G_PP)
    return m_pi, G_PP, f_pi


def _f_pi_from_pcac(m_q, m_pi, G_PP):
    """f_pi*a (bare/unrenormalised) from the PCAC / pseudoscalar-density route:
        f_pi = (2 m_q / m_pi^2) * sqrt(G_PP),
    with m_q the bare Wilson quark mass (input mass shifted by the critical mass) and
    G_PP = |<0|Pbar|pi>|^2 the pseudoscalar-density overlap. This is the unrenormalised number;
    the renormalised f_pi needs Z_A, Z_P and m_crit (all OWED). Convention note: this returns
    f_pi in the f_pi ~ 92 MeV (chiral, F_pi) normalisation; multiply by sqrt(2) for the 130 MeV
    (f_pi = sqrt(2) F_pi) convention."""
    if not (np.isfinite(m_pi) and np.isfinite(G_PP)) or m_pi <= 0 or G_PP < 0:
        return float("nan")
    return float((2.0 * m_q / (m_pi ** 2)) * np.sqrt(G_PP))


def decay_constant(c_pp, m_q, T, tmin=None, tmax=None):
    r"""The pion decay constant f_pi*a (and the substrate scale v) from the zero-momentum
    pseudoscalar correlator C_PP(t) (target L2: test sigma = 2 pi v^2).

    `c_pp` is the (n,3) array of rows [cfg, t, C_PP(t)] (measure_decay_constant output, tagged per
    config), C_PP(t) = sum_x Tr[gamma5 S(x,t;0) gamma5 S(x,t;0)^dag] = sum_x Tr[S^dag S]. The config
    average is fit to the two-sided cosh
        C_PP(t) = (G_PP / (2 m_pi)) * (e^{-m_pi t} + e^{-m_pi (T-t)})
    to get m_pi*a and the pseudoscalar overlap G_PP = |<0|Pbar|pi>|^2, then f_pi*a via the PCAC /
    pseudoscalar-density relation
        f_pi = (2 m_q / m_pi^2) * sqrt(G_PP)     (bare/unrenormalised; m_q = bare Wilson quark mass).
    Errors are delete-1 jackknife over configs.

    Substrate scale v: in the framework's NJL/chiral-soliton identification v = f_pi (this is an
    ASSUMPTION of the framework, stated not hidden), so v*a = f_pi*a here; the alternative
    v^2 ~ |<qbar q>| is not computed here.

    HONEST: single (heavy) sea mass; UNRENORMALISED (Z_A, Z_P, m_crit all owed); the f_pi ~ 92 MeV
    (F_pi) normalisation is used (multiply by sqrt(2) for the 130 MeV convention). Returns m_pi*a,
    f_pi*a (= v*a), G_PP, jackknife errors, the fit window, and the per-config table."""
    c_pp = np.asarray(c_pp, dtype=float)
    cfg = c_pp[:, 0].astype(int)
    tt = c_pp[:, 1].astype(int)
    C = c_pp[:, 2]
    ts = np.arange(tt.max() + 1)
    configs = np.unique(cfg)
    n = len(configs)

    Cbar = np.array([C[tt == s].mean() for s in ts])
    m_pi, G_PP, f_pi = _cosh_fit_pp(Cbar, m_q, T, tmin, tmax)

    jm, jg, jf = [], [], []
    for cc in configs:                                   # delete-1 jackknife over configs
        keep = cfg != cc
        Cj = np.array([C[keep & (tt == s)].mean() for s in ts])
        mj, gj, fj = _cosh_fit_pp(Cj, m_q, T, tmin, tmax)
        if np.isfinite(mj):
            jm.append(mj)
        if np.isfinite(gj):
            jg.append(gj)
        if np.isfinite(fj):
            jf.append(fj)

    def jk_err(vals):
        vals = np.asarray(vals, dtype=float)
        if vals.size < 2:
            return float("nan")
        return float(np.sqrt((vals.size - 1) / vals.size * np.sum((vals - vals.mean()) ** 2)))

    lo = tmin if tmin is not None else max(1, T // 8)
    hi = tmax if tmax is not None else max(lo + 1, T // 2 - 1)
    return {"m_pi": m_pi, "m_pi_err": jk_err(jm), "G_PP": G_PP, "G_PP_err": jk_err(jg),
            "f_pi": f_pi, "f_pi_err": jk_err(jf), "v": f_pi, "v_err": jk_err(jf),
            "m_q": float(m_q), "tmin": int(lo), "tmax": int(hi), "n_cfg": int(n),
            "convention": "f_pi ~ 92 MeV (F_pi); x sqrt(2) for the 130 MeV convention",
            "renorm_note": "UNRENORMALISED: Z_A, Z_P, m_crit owed; v=f_pi is the framework NJL assumption"}


def sigma_2piv2_check(sqrt_sigma_a, f_pi_a, sqrt_sigma_a_err=None, f_pi_a_err=None, tol=0.5):
    r"""Test the framework's sharp prediction sigma = 2 pi v^2 (target L2): one substrate sets BOTH
    the mass scale and the confinement string tension. With v identified as f_pi (the NJL/chiral-
    soliton assumption), the lattice numbers are sqrt(sigma)*a (static potential) and f_pi*a (this
    measurement). The dimensionless test is the ratio

        ratio = sigma / (2 pi v^2) = (sqrt(sigma)*a)^2 / (2 pi (f_pi*a)^2),

    which the prediction sets to 1 (the lattice spacing a cancels). `tol` is the fractional band
    around 1 within which we call it consistent (default 0.5: a factor-2 band, appropriate for a
    single heavy sea mass with unrenormalised f_pi). Returns the ratio, its error (if inputs carry
    errors), consistency flag, and a verdict STRING.

    HONEST: single heavy sea mass; f_pi*a is UNRENORMALISED (Z_A, Z_P, m_crit owed) and the v = f_pi
    identification is the framework's NJL assumption, not a derived equality -- so a ratio off 1 by a
    Z-factor (~O(1)) is expected and NOT yet a falsification; only a gross (order-of-magnitude)
    mismatch would be."""
    sa = float(sqrt_sigma_a)
    fa = float(f_pi_a)
    sigma = sa ** 2
    twopi_v2 = 2.0 * np.pi * fa ** 2
    ratio = float(sigma / twopi_v2) if twopi_v2 != 0 else float("nan")

    # error propagation: ratio = sigma_a^2 / (2 pi f_a^2); d ln ratio = 2 dsa/sa - 2 dfa/fa
    ratio_err = float("nan")
    if sqrt_sigma_a_err is not None and f_pi_a_err is not None and np.isfinite(ratio):
        rel = np.sqrt((2.0 * sqrt_sigma_a_err / sa) ** 2 + (2.0 * f_pi_a_err / fa) ** 2)
        ratio_err = float(abs(ratio) * rel)

    consistent = bool(np.isfinite(ratio) and abs(ratio - 1.0) <= tol)
    # is it consistent within the propagated error band (a looser, statistics-aware test)?
    within_err = bool(np.isfinite(ratio_err) and abs(ratio - 1.0) <= ratio_err)

    f_pi_pred_a = float(np.sqrt(sigma / (2.0 * np.pi))) if sigma >= 0 else float("nan")

    if not np.isfinite(ratio):
        verdict = ("sigma = 2 pi v^2 UNTESTABLE here: f_pi*a did not come out finite (no pion plateau "
                   "/ amplitude). Get a clean C_PP cosh fit first.")
    elif consistent:
        verdict = (f"sigma = 2 pi v^2 CONSISTENT: sigma/(2 pi f_pi^2) = {ratio:.2f}"
                   + (f" +/- {ratio_err:.2f}" if np.isfinite(ratio_err) else "")
                   + f" (within the factor-{1+tol:.0f} band of 1). One substrate plausibly sets both the "
                   f"mass scale and the string tension: from sqrt(sigma)*a={sa:.3f} the prediction is "
                   f"f_pi*a={f_pi_pred_a:.3f} vs measured {fa:.3f}. CAVEATS: single heavy sea mass; f_pi*a "
                   f"is UNRENORMALISED (Z_A, Z_P, m_crit all owed, each an O(1) shift); v=f_pi is the "
                   f"framework's NJL assumption, not a derived equality -- so this is an encouraging "
                   f"order-of-magnitude success, not yet a confirmed equality.")
    elif within_err:
        verdict = (f"sigma = 2 pi v^2 consistent WITHIN ERRORS: sigma/(2 pi f_pi^2) = {ratio:.2f} +/- "
                   f"{ratio_err:.2f} brackets 1, though the central value is off the factor-{1+tol:.0f} "
                   f"band. Same caveats (single heavy mass; unrenormalised; v=f_pi assumed). Tighten the "
                   f"statistics / add the Z-factors before reading the central offset as physics.")
    else:
        direction = "LARGER" if ratio > 1.0 else "SMALLER"
        verdict = (f"sigma = 2 pi v^2 STRAINED: sigma/(2 pi f_pi^2) = {ratio:.2f}"
                   + (f" +/- {ratio_err:.2f}" if np.isfinite(ratio_err) else "")
                   + f" -- sigma is {direction} than 2 pi f_pi^2 by more than the factor-{1+tol:.0f} band "
                   f"(predicted f_pi*a={f_pi_pred_a:.3f} vs measured {fa:.3f}). NOT yet a falsification: "
                   f"the f_pi*a here is UNRENORMALISED (a Z_A/Z_P/m_crit factor of O(1) is owed) at a "
                   f"single heavy sea mass, and v=f_pi is the framework's NJL assumption. Add the "
                   f"renormalisation and a chiral extrapolation before reading this as a verdict against "
                   f"the prediction; only an order-of-magnitude mismatch would falsify.")

    return {"ratio": ratio, "ratio_err": ratio_err, "consistent": consistent,
            "within_err": within_err, "sqrt_sigma_a": sa, "f_pi_a": fa, "sigma_a2": sigma,
            "twopi_v2_a2": twopi_v2, "f_pi_pred_a": f_pi_pred_a, "tol": float(tol),
            "verdict": verdict}


# ---------------------------------------------------------------------------
# 5b. The electroweak S parameter proxy from the isovector V-A correlators (the L3 control)
# ---------------------------------------------------------------------------
def _cosh_fit_ga(Cbar, T, tmin=None, tmax=None):
    """Two-sided cosh fit of a zero-momentum meson correlator C(t) to
        C(t) = (G / (2 M)) * (e^{-M t} + e^{-M (T-t)}),
    returning (M*a, G) -- the mass from the periodic acosh effective-mass plateau and the amplitude
    G from the least-squares overlap of the fixed-shape cosh, exactly as _cosh_fit_pp but WITHOUT the
    PCAC f_pi step (the vector/axial channels carry their own decay-constant convention, handled in
    s_parameter_proxy). G is the SAME amplitude convention as _cosh_fit_pp: C(t)=(G/(2M))*shape(t)."""
    C = np.asarray(Cbar, dtype=float)
    if tmin is None:
        tmin = max(1, T // 8)
    if tmax is None:
        tmax = max(tmin + 1, T // 2 - 1)
    tt = np.arange(len(C))
    meff = np.full(len(C), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        for t in range(1, len(C) - 1):
            num = C[t - 1] + C[t + 1]
            den = 2.0 * C[t]
            if den != 0 and num / den >= 1.0:
                meff[t] = np.arccosh(num / den)
    idx = np.arange(tmin, min(tmax, len(C) - 1))
    idx = idx[np.isfinite(meff[idx])]
    if len(idx) < 1:
        return float("nan"), float("nan")
    M = float(np.mean(meff[idx]))
    shape = np.exp(-M * tt) + np.exp(-M * (T - tt))
    A = (shape[idx] / (2.0 * M))
    denom = float(np.dot(A, A))
    G = float(np.dot(A, C[idx]) / denom) if denom > 0 else float("nan")
    return M, G


def _f_from_G(G, M, npol=3):
    """The vector/axial decay constant F (per polarisation) from the cosh amplitude G of the
    zero-momentum, spatial-polarisation-SUMMED correlator. With <0|V_i|V> = F M eps_i and the
    correlator summed over npol=3 spatial polarisations, the cosh amplitude is G = npol * F^2 * M^2,
    so F^2 = G / (npol * M^2) and F = sqrt(G/(npol*M^2)). Returns F (NaN if G<0 or M<=0)."""
    if not (np.isfinite(G) and np.isfinite(M)) or M <= 0 or G < 0:
        return float("nan")
    return float(np.sqrt(G / (npol * M ** 2)))


def s_parameter_proxy(c_v, c_a, T, tmin=None, tmax=None, npol=3):
    r"""A FIRST lattice estimate of the electroweak S parameter (Peskin-Takeuchi) from the isovector
    VECTOR minus AXIAL-VECTOR correlators (measure_s_parameter), the L3 control point (Eternal Dawn,
    Ch. "The Forces from the Field": does the electroweak-breaking sector clear S < 0.1, or sit in the
    QCD-like graveyard at S ~ 0.25?).

    Inputs (measure_s_parameter rows, config-tagged):
      c_v : (n,3) [cfg, t, C_V(t)]   -- the connected isovector VECTOR correlator (gamma_i, rho channel)
      c_a : (n,3) [cfg, t, C_A(t)]   -- the connected isovector AXIAL correlator (gamma_i g5, a1 channel)
    both summed over the npol=3 spatial polarisations i=1,2,3 (the disconnected pieces cancel for the
    isovector current).

    Method. The config-averaged C_V and C_A are each cosh-fit (the SAME machinery as decay_constant) for
    the rho mass M_V, the a1 mass M_A, and the cosh amplitudes G_V, G_A; the decay constants follow from
    F^2 = G/(npol M^2) (so F^2/M^2 = G/(npol M^4)). The lowest-resonance (vector-meson-dominance /
    Weinberg-sum-rule-saturated) S-proxy is

        S = 4 pi ( F_V^2 / M_V^2  -  F_A^2 / M_A^2 ).

    Also returned: the DIRECT walking diagnostic M_A/M_V (-> 1 if conformal, ~1.6 QCD-like) and the
    integrated chiral order parameter sum_t (C_V(t) - C_A(t)) (nonzero iff chiral symmetry is broken).
    Errors are delete-1 jackknife over configs.

    HONEST (mirror the bar of decay_constant / condensate_3pt): this is the LOWEST-RESONANCE saturation
    (the rho and a1 only, not the full spectral integral); the currents are UNRENORMALISED (Z_V, Z_A are
    owed, each an O(1) factor); a SINGLE heavy sea mass (no chiral extrapolation); and the FUNDAMENTAL
    rep is QCD-like, NOT walking, so this is EXPECTED to land near S ~ 0.25. That characterises the
    sector and TESTS the assumption -- it does NOT validate S < 0.1, which needs the near-conformal /
    propagating-torsion regime (target L5). Returns M_V, M_A, M_A/M_V, F_V, F_A, S, integrated V-A, the
    jackknife errors, the fit window, and a verdict STRING."""
    c_v = np.asarray(c_v, dtype=float)
    c_a = np.asarray(c_a, dtype=float)

    def _avg(raw):
        cfg = raw[:, 0].astype(int)
        tt = raw[:, 1].astype(int)
        C = raw[:, 2]
        ts = np.arange(tt.max() + 1)
        Cbar = np.array([C[tt == s].mean() for s in ts])
        return cfg, tt, C, ts, Cbar

    cfgV, ttV, CV, tsV, CVbar = _avg(c_v)
    cfgA, ttA, CA, tsA, CAbar = _avg(c_a)
    configs = np.unique(np.concatenate([cfgV, cfgA]))
    n = len(configs)

    def _point(CVb, CAb):
        M_V, G_V = _cosh_fit_ga(CVb, T, tmin, tmax)
        M_A, G_A = _cosh_fit_ga(CAb, T, tmin, tmax)
        F_V = _f_from_G(G_V, M_V, npol)
        F_A = _f_from_G(G_A, M_A, npol)
        S = float("nan")
        if all(np.isfinite(x) for x in (F_V, M_V, F_A, M_A)) and M_V > 0 and M_A > 0:
            S = 4.0 * np.pi * (F_V ** 2 / M_V ** 2 - F_A ** 2 / M_A ** 2)
        ratio = float(M_A / M_V) if (np.isfinite(M_V) and M_V > 0) else float("nan")
        # integrated chiral order parameter sum_t (C_V - C_A) over the common time extent
        ncom = min(len(CVb), len(CAb))
        va_int = float(np.sum(CVb[:ncom] - CAb[:ncom]))
        return {"M_V": M_V, "M_A": M_A, "ratio": ratio, "F_V": F_V, "F_A": F_A,
                "S": S, "va_int": va_int}

    full = _point(CVbar, CAbar)

    jk = {k: [] for k in ("M_V", "M_A", "ratio", "F_V", "F_A", "S", "va_int")}
    for cc in configs:                                   # delete-1 jackknife over configs
        kv = cfgV != cc
        ka = cfgA != cc
        CVj = np.array([CV[kv & (ttV == s)].mean() for s in tsV])
        CAj = np.array([CA[ka & (ttA == s)].mean() for s in tsA])
        pj = _point(CVj, CAj)
        for k in jk:
            if np.isfinite(pj[k]):
                jk[k].append(pj[k])

    def jk_err(vals):
        vals = np.asarray(vals, dtype=float)
        if vals.size < 2:
            return float("nan")
        return float(np.sqrt((vals.size - 1) / vals.size * np.sum((vals - vals.mean()) ** 2)))

    lo = tmin if tmin is not None else max(1, T // 8)
    hi = tmax if tmax is not None else max(lo + 1, T // 2 - 1)

    S, S_err = full["S"], jk_err(jk["S"])
    M_V, M_A, ratio = full["M_V"], full["M_A"], full["ratio"]
    va_int = full["va_int"]
    Q_QCD, S_TARGET = 0.25, 0.10                          # the QCD graveyard value and the walking target

    if not np.isfinite(S):
        verdict = ("S-proxy UNTESTABLE here: a clean cosh fit for M_V and/or M_A did not come out (no "
                   "vector / axial plateau or amplitude). Get clean C_V and C_A cosh fits first.")
    else:
        # the SIGN check first: a positive integrated V-A (chiral symmetry broken) and S>0 is the
        # QCD-like expectation; M_A/M_V well above 1 (~1.6) confirms the a1 sits above the rho (not
        # near-conformal). State plainly that the fundamental rep is EXPECTED here.
        if S > 0:
            near = abs(S - Q_QCD) <= 0.5 * Q_QCD          # within a factor ~1.5 of the QCD value
            place = ("near the QCD value ~0.25 (the EXPECTED QCD-like result)" if near
                     else (f"ABOVE the QCD value ~0.25" if S > Q_QCD
                           else f"BELOW the QCD value ~0.25 but still positive"))
            verdict = (
                f"S-proxy = {S:.3f}" + (f" +/- {S_err:.3f}" if np.isfinite(S_err) else "")
                + f", {place}. The walking diagnostic M_A/M_V = {ratio:.2f} "
                + ("(QCD-like, a1 clearly above rho)" if np.isfinite(ratio) and ratio > 1.3
                   else "(close to 1 -- check the fit, this would be the conformal signature)")
                + f"; the integrated V-A = {va_int:.3e} "
                + ("(positive: chiral symmetry broken, as expected)" if va_int > 0
                   else "(non-positive -- recheck the channels)")
                + f". HONEST: this is the lowest-resonance (VMD/Weinberg-saturated) proxy, currents "
                + f"UNRENORMALISED (Z_V, Z_A owed), a SINGLE heavy sea mass, and the FUNDAMENTAL rep is "
                + f"QCD-like, NOT walking -- so landing near {Q_QCD:.2f} CHARACTERISES the sector and "
                + f"TESTS the load-bearing assumption; it does NOT validate S < {S_TARGET:.2f}, which "
                + f"needs the near-conformal / propagating-torsion regime (target L5)."
            )
        else:
            verdict = (
                f"S-proxy = {S:.3f}" + (f" +/- {S_err:.3f}" if np.isfinite(S_err) else "")
                + f" is NEGATIVE -- below the naive QCD-like expectation. With M_A/M_V = {ratio:.2f} and "
                + f"integrated V-A = {va_int:.3e}, a negative lowest-resonance S can come from an "
                + f"under-resolved a1 (heavy, noisy axial channel) or the missing higher-resonance / "
                + f"contact saturation, NOT from genuine walking. UNRENORMALISED (Z_V, Z_A owed), single "
                + f"heavy sea mass; this does NOT validate S < {S_TARGET:.2f} (that needs the near-"
                + f"conformal L5 regime). Tighten the axial fit / add resonances before reading the sign."
            )

    return {"M_V": M_V, "M_V_err": jk_err(jk["M_V"]), "M_A": M_A, "M_A_err": jk_err(jk["M_A"]),
            "ratio": ratio, "ratio_err": jk_err(jk["ratio"]),
            "F_V": full["F_V"], "F_V_err": jk_err(jk["F_V"]),
            "F_A": full["F_A"], "F_A_err": jk_err(jk["F_A"]),
            "S": S, "S_err": S_err, "va_int": va_int, "va_int_err": jk_err(jk["va_int"]),
            "tmin": int(lo), "tmax": int(hi), "n_cfg": int(n), "npol": int(npol),
            "S_qcd_ref": Q_QCD, "S_target": S_TARGET,
            "renorm_note": "UNRENORMALISED: Z_V, Z_A owed; lowest-resonance saturation; single heavy sea "
                           "mass; fundamental rep EXPECTED QCD-like (~0.25), NOT walking -- does not "
                           "validate S<0.1 (needs near-conformal/L5).",
            "verdict": verdict}


def bag_chiral_trend(points, window=(0.43, 0.70), r0_over_a=3.166):
    """Chiral extrapolation of the bag sharpness s_T=R0/r0 vs the pion mass squared (Eternal Dawn,
    Ch. 11 Generations -- the real test, not a single heavy point).

    `points` is a list of (m_pi2, s_T, s_T_err) across valence masses (m_pi2 in lattice units (m_pi a)^2;
    s_T = R0/r0 dimensionless). A heavy quark gives a SMALL bag (small s_T); the question is whether,
    as m_pi2 -> 0, s_T RISES into the lever's productive window [0.43,0.70] r0 (where the consecutive
    ladder reproduces the observed lepton span ~3477). We fit s_T = c0 + c1 m_pi2 (weighted by 1/err^2)
    and read the chiral intercept c0 = s_T(m_pi=0).

    Returns the chiral s_T and error, the slope, whether the trend RISES toward chiral (c1<0), whether
    the chiral bag lands IN the window, the implied span there (only quoted if in-window), and a
    verdict. HONEST: this is the VALENCE (partially-quenched) chiral limit at fixed sea; the unitary
    limit and the absolute lever normalisation are the further steps."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts.reshape(1, -1)
    x = pts[:, 0]
    y = pts[:, 1]
    e = pts[:, 2] if pts.shape[1] > 2 else np.full_like(y, np.nan)
    lo, hi = window

    if len(x) < 2:
        return {"chiral_s_T": float(y[0]), "chiral_s_T_err": float("nan"), "slope": float("nan"),
                "rising": None, "in_window": bool(lo <= y[0] <= hi), "span": float("nan"),
                "points": pts, "fit": (float(y[0]), 0.0),
                "verdict": "single mass only -- need >=2 valence points to extrapolate the chiral trend."}

    w = 1.0 / np.clip(e, 1e-9, None) ** 2
    if not np.all(np.isfinite(w)):
        w = np.ones_like(y)

    def _wlinfit(xx, yy, ww):
        Ak = np.vstack([np.ones_like(xx), xx]).T
        AtW = Ak.T * ww
        covk = np.linalg.inv(AtW @ Ak)
        ck = covk @ (AtW @ yy)
        return float(ck[0]), float(np.sqrt(covk[0, 0])), float(ck[1])

    c0, c0_err, c1 = _wlinfit(x, y, w)                  # all-point fit (conservative)
    rising = bool(c1 < 0)                               # s_T grows as m_pi2 falls toward chiral
    in_window = bool(lo <= c0 <= hi)
    lightest_s_T = float(y[np.argmin(x)])               # s_T at the lightest (most chiral) point
    any_in = bool(np.any((y >= lo) & (y <= hi)))        # is any MEASURED point already in the window?
    n_in = int(np.sum((y >= lo) & (y <= hi)))
    accel = bool(len(x) >= 3 and rising and lightest_s_T > c0 + c0_err)

    # CHIRAL extrapolation should use the LIGHT points (heavy masses curve out of the chiral regime).
    # Report the fit as the heaviest points are dropped progressively, lightest-first -- transparent,
    # not cherry-picked -- and take the light-half fit as the chiral estimate.
    order = np.argsort(x)
    xs, ys, ws = x[order], y[order], w[order]
    stability = []                                      # (n_points_used, c0, c0_err), lightest-first
    for k in range(len(xs), 2, -1):
        ck0, cke, _ = _wlinfit(xs[:k], ys[:k], ws[:k])
        stability.append((k, ck0, cke))
    k_light = max(3, len(xs) // 2)
    c0_light, c0_light_err = next((c, ce) for (kk, c, ce) in stability if kk == k_light)
    in_window_light = bool(lo <= c0_light <= hi)

    span = float("nan")
    s_for_span = (c0_light if in_window_light else (lightest_s_T if any_in else c0))
    if in_window_light or any_in or in_window:
        from . import fermion_masses as fm
        lad = fm._ascending_ladder(float(np.clip(s_for_span, 0.30, 1.5)), well=fm.WELL)
        span = float(lad[-1] / lad[0])

    if in_window or in_window_light:
        verdict = (f"the bag GROWS into the productive window. Chiral extrapolation: all-mass "
                   f"{c0:.2f}+/-{c0_err:.2f} r0 (a curvature-contaminated LOWER bound), light-mass "
                   f"(lightest {k_light}) {c0_light:.2f}+/-{c0_light_err:.2f} r0 IN [{lo},{hi}]; "
                   f"{n_in} measured point(s) already in -> lever span ~{span:.0f} vs observed 3477. "
                   f"CONSISTENT WITH DERIVED (valence chiral limit; confirm with the 3-pt condensate).")
    elif rising and any_in:
        verdict = (f"ENCOURAGING: light-mass chiral extrapolation s_T={c0_light:.2f} r0 (all-mass "
                   f"{c0:.2f}, a lower bound), {n_in} measured point(s) already IN the window "
                   f"(lightest s_T={lightest_s_T:.2f}, span ~{span:.0f})"
                   + (" and the rise ACCELERATES toward chiral" if accel else "")
                   + ". Firm with more light-mass stats / the 3-pt condensate.")
    elif rising:
        verdict = (f"the bag RISES toward chiral but the light-mass extrapolation s_T={c0_light:.2f} r0 "
                   f"is still BELOW the window [{lo},{hi}] -- the one-body bag under-shoots; push to "
                   f"lighter mass and measure the 3-pt condensate <N|qbar q(r)|N>.")
    elif not rising:
        verdict = (f"the bag does NOT grow toward chiral (slope {c1:+.2f}); chiral s_T = {c0:.2f} r0. "
                   f"The one-body proxy under-shoots -- the hierarchy needs the 3-pt condensate.")
    else:
        verdict = (f"chiral s_T = {c0:.2f}+/-{c0_err:.2f} r0 ABOVE the window [{lo},{hi}] -- bag over-"
                   f"broad; recheck the scale and window.")

    return {"chiral_s_T": float(c0), "chiral_s_T_err": c0_err, "slope": float(c1), "rising": rising,
            "in_window": in_window, "span": span, "points": pts, "fit": (float(c0), float(c1)),
            "lightest_s_T": lightest_s_T, "any_in_window": any_in, "n_in_window": n_in, "accelerating": accel,
            "chiral_s_T_light": float(c0_light), "chiral_s_T_light_err": float(c0_light_err),
            "k_light": int(k_light), "in_window_light": in_window_light, "stability": stability,
            "verdict": verdict}


def _twostate_profile(p3, sr, t_snk, r0_over_a, taumin=1):
    """Ground-state condensate profile from a TWO-STATE fit in the insertion time tau (the right tool
    for the SHAPE: it removes the time-direction excited-state contamination of a POINT source without
    the spatial broadening that source smearing imposes). At fixed sink, the 3-pt ratio has the form
        G3(r,tau) = M(r) + A(r) * [ e^{-dE*tau} + e^{-dE*(t_snk-tau)} ],
    M(r) the ground-state matrix element. We fix the gap dE from the integrated 3-pt SR(tau) (best
    signal/noise), then a per-shell LINEAR fit gives M(r). Returns (r, rho_gs, R0, s_T, dE) or None."""
    try:
        from scipy.optimize import curve_fit
    except Exception:
        return None
    p3 = np.asarray(p3, dtype=float)
    sr = np.asarray(sr, dtype=float)
    taus = np.arange(taumin, t_snk - taumin + 1)
    if len(taus) < 4:
        return None
    SRbar = np.array([sr[sr[:, 1].astype(int) == t, 2].mean() if np.any(sr[:, 1].astype(int) == t)
                      else np.nan for t in taus])
    if not np.all(np.isfinite(SRbar)):
        return None
    g = lambda tau, dE: np.exp(-dE * tau) + np.exp(-dE * (t_snk - tau))
    try:                                                  # gap dE from the integrated 3-pt
        popt, _ = curve_fit(lambda tau, a, b, dE: a + b * g(tau, dE), taus, SRbar,
                            p0=[SRbar.mean(), float(np.ptp(SRbar)) or 1.0, 0.5],
                            bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 3.0]), maxfev=20000)
        dE = float(popt[2])
    except Exception:
        return None
    gvec = g(taus.astype(float), dE)
    # per-shell linear fit M(r) + A(r) gvec, aggregated over configs
    r2s = np.unique(p3[:, 1]).astype(int)
    r, M = [], []
    for b in r2s:
        rows = p3[p3[:, 1].astype(int) == b]
        y = np.array([rows[rows[:, 2].astype(int) == t, 3].sum() /
                      max(rows[rows[:, 2].astype(int) == t, 4].sum(), 1.0)
                      for t in taus])
        if not np.all(np.isfinite(y)):
            continue
        A = np.vstack([np.ones_like(gvec), gvec]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r.append(np.sqrt(b)); M.append(coef[0])           # coef[0] = ground-state M(r)
    r, M = np.asarray(r, float), np.asarray(M, float)
    order = np.argsort(r); r, M = r[order], M[order]
    if len(M) < 3 or M[0] == 0:
        return None
    rho_gs = M / M[0]
    R0 = _half_density_radius(r, rho_gs)
    return {"r": r, "rho_gs": rho_gs, "R0": float(R0), "s_T": float(R0 / r0_over_a), "dE": dE}


def condensate_3pt(p3, c2, sr, chk, T, t_snk, tau_window=None, r0_over_a=3.166, selfcheck_tol=2e-2):
    """The connected nucleon scalar 3-point <N|qbar q(r)|N> -- the genuine in-medium condensate bag
    (measure_condensate_3pt), the definitive s_T (Eternal Dawn, Ch. 11 Generations).

    Inputs (parsed measure_condensate_3pt rows, config-tagged):
      p3   : (n,5) [cfg, r2, tau, G3_sum, count]   -- the shell-summed 3-point density
      c2   : (n,3) [cfg, t, C_N(t)]                -- the 2-point (for the sink normalisation)
      sr   : (n,3) [cfg, tau, sum_x G3(tau)]       -- the sum-rule numerator
      chk  : (n,3) [cfg, recon_2pt, C_N(t_snk)]    -- the sequential-source SELF-CHECK pair

    THE GATE: the self-check. C_N is linear in the spectator propagator, so the source reconstruction
    MUST equal C_N(t_snk). If |recon - C_N|/|C_N| > selfcheck_tol the sequential SOURCE is wrong and
    nothing downstream is trusted -- the verdict says so first. If it passes: the sum rule
    g_S = <sum_x G3(tau)/C_N(t_snk)> plateaus in tau (the nucleon scalar charge), and the tau-plateau
    profile G3(r) is the condensate bag -- its half-density radius gives s_T = R0/r0 DIRECTLY (no one-
    body dictionary). Returns the self-check status, g_S, the profile, R0, s_T (with jackknife), and a
    verdict leading with the gate."""
    p3 = np.asarray(p3, dtype=float)
    c2 = np.asarray(c2, dtype=float)
    sr = np.asarray(sr, dtype=float)
    chk = np.asarray(chk, dtype=float)

    # --- the self-check (config-averaged reconstruction vs C_N(t_snk)) ---
    recon = float(chk[:, 1].mean())
    cn_snk = float(chk[:, 2].mean())
    sc_resid = abs(recon - cn_snk) / (abs(cn_snk) + 1e-300)
    sc_ok = bool(sc_resid <= selfcheck_tol)

    # --- sink quality: the nucleon 2pt must be on the plateau (forward state dominant) at t_snk.
    # The contraction's OVERALL SIGN is a convention (C_N can be uniformly negative), so judge by the
    # SIGN-NORMALISED correlator: positive & decaying on the forward side, turning over near T/2 where
    # the backward-propagating wrong-parity state takes over (the real node). ---
    tt2 = c2[:, 1].astype(int)
    c2bar = np.array([c2[tt2 == t, 2].mean() if np.any(tt2 == t) else np.nan
                      for t in range(int(tt2.max()) + 1)])
    ref = c2bar[2:7][np.isfinite(c2bar[2:7])]            # signal region, past the contact term
    sgn = -1.0 if (len(ref) and np.sum(ref) < 0) else 1.0
    sc2 = sgn * c2bar                                    # sign-normalised: forward side > 0, decaying
    node_t = next((t for t in range(2, len(sc2)) if not (sc2[t] > 0)), len(sc2))
    sink_ok = bool(sgn * cn_snk > 0 and t_snk < node_t)

    # --- the sum rule g_S(tau) = sum_x G3(tau) / C_N(t_snk), plateau between source and sink ---
    taus = np.arange(int(sr[:, 1].max()) + 1)
    SRbar = np.array([sr[sr[:, 1].astype(int) == t, 2].mean() if np.any(sr[:, 1].astype(int) == t)
                      else np.nan for t in taus])
    gS_tau = SRbar / (cn_snk + 1e-300)
    if tau_window is None:
        tau_window = (max(1, t_snk // 4), max(2, t_snk - t_snk // 4))
    lo, hi = tau_window
    win = np.arange(lo, min(hi, len(gS_tau)))
    win = win[np.isfinite(gS_tau[win])]
    g_S = float(np.mean(gS_tau[win])) if len(win) else float("nan")

    # --- the tau-plateau condensate profile G3(r) and its half-density radius ---
    r, rho3 = _bag_rho(p3, tau_window)                  # reuse: aggregate shells over the tau window
    R0 = _half_density_radius(r, rho3)
    s_T = R0 / r0_over_a

    # ground-state shape via the two-state fit in tau (the right tool -- point source, no smearing)
    gs = _twostate_profile(p3, sr, t_snk, r0_over_a)
    gs_s_T = gs["s_T"] if gs else float("nan")
    gs_R0 = gs["R0"] if gs else float("nan")

    cfgs = np.unique(p3[:, 0]).astype(int)
    n = len(cfgs)
    js = []
    for c in cfgs:
        rj, rhoj = _bag_rho(p3[p3[:, 0].astype(int) != c], tau_window)
        js.append(_half_density_radius(rj, rhoj) / r0_over_a)
    js = np.array(js)
    s_T_err = float(np.sqrt((n - 1) / n * np.sum((js - js.mean()) ** 2))) if n > 1 else float("nan")

    lo_p, hi_p = 0.43, 0.70
    where = ("IN" if lo_p <= s_T <= hi_p else ("below" if s_T < lo_p else "above"))
    # the two-state fit is only TRUSTWORTHY when its R0 is comparable to the (resolved) plateau R0; if
    # it swings far from it the tau-fit is ill-conditioned (it over/under-corrected). At a large t_snk
    # the plateau is itself ground-state-dominated, so prefer it and flag the two-state as unstable.
    gs_reliable = bool(gs is not None and R0 > 1.2 and 0.6 < (gs_R0 / max(R0, 1e-9)) < 2.2)
    if not sc_ok:
        verdict = (f"SELF-CHECK FAILED: reconstruction {recon:.4g} vs C_N(t_snk) {cn_snk:.4g} "
                   f"(resid {sc_resid:.1%}) -- the sequential source is wrong; fix the sigma_seq "
                   f"sign/transpose before trusting g_S or the profile.")
    elif not sink_ok:
        verdict = (f"self-check PASSED (the contraction is correct), BUT the SINK IS PAST THE PLATEAU: "
                   f"the sign-normalised 2-pt turns over at t={node_t} (backward-propagating wrong-"
                   f"parity state), and t_snk={t_snk} is at/beyond it. The s_T={s_T:.3f} here is "
                   f"backward-contaminated; rerun with SINKT < {node_t}.")
    elif R0 > 1.2 and not gs_reliable:
        # RESOLVED bag, large t_snk -> the plateau IS the ground state; the two-state fit is unstable
        twostate = f" (two-state fit unstable here: R0={gs_R0:.1f}a, ignore it)" if gs is not None else ""
        verdict = (f"self-check PASSED, sink on the plateau (t_snk={t_snk} < node {node_t}); g_S = "
                   f"{abs(g_S):.3f}. RESOLVED condensate bag s_T = R0/r0 = {s_T:.3f}+/-{s_T_err:.3f} r0 "
                   f"({where} the window [{lo_p},{hi_p}]) -- the genuine three-body number, R0={R0:.2f}a "
                   f"off the cutoff, plateau ground-state-dominated at this large t_snk{twostate}.")
    elif gs is not None:
        gw = ("IN" if lo_p <= gs_s_T <= hi_p else ("below" if gs_s_T < lo_p else "above"))
        verdict = (f"self-check PASSED, sink on the plateau (t_snk={t_snk} < node {node_t}); g_S = "
                   f"{abs(g_S):.3f}. GROUND-STATE bag (two-state fit in tau, dE={gs['dE']:.2f}): "
                   f"s_T = R0/r0 = {gs_s_T:.3f} r0 ({gw} the window [{lo_p},{hi_p}]). "
                   f"[the tau-plateau average {s_T:.3f} is excited-contaminated -- raise t_snk so the "
                   f"plateau is ground-state-dominated.]")
    else:
        verdict = (f"self-check PASSED (resid {sc_resid:.1%}), sink on the plateau (t_snk={t_snk} < "
                   f"node {node_t}); scalar charge g_S = {abs(g_S):.3f}. tau-plateau bag s_T = R0/r0 = "
                   f"{s_T:.3f}+/-{s_T_err:.3f} r0 ({where} the window [{lo_p},{hi_p}]). NOTE: point-source "
                   f"3-pt, plateau average -- excited-contaminated; the two-state fit needs more tau "
                   f"points (bigger T) to isolate the ground state.")

    # the best single number: the resolved plateau when the two-state fit is unstable, else the gs fit
    s_T_best = s_T if (R0 > 1.2 and not gs_reliable) else (gs_s_T if gs is not None else s_T)
    return {"self_check_ok": sc_ok, "self_check_resid": sc_resid, "recon": recon, "cn_snk": cn_snk,
            "sink_ok": sink_ok, "node_t": int(node_t), "g_S": g_S, "g_S_tau": gS_tau, "r": r,
            "rho3": rho3, "R0": R0, "s_T": s_T, "s_T_err": s_T_err, "n_cfg": n,
            "gs_s_T": gs_s_T, "gs_R0": gs_R0, "gs": gs, "gs_reliable": gs_reliable, "s_T_best": s_T_best,
            "tau_window": tau_window, "verdict": verdict}


def _gevp_lambdas(Cmat, t0, eps=1e-6):
    """Robust generalized eigenvalues lambda_n(t) of C(t) v = lambda C(t0) v, sorted DESCENDING.
    Rather than require C(t0) positive-definite (it rarely is, once a point operator is mixed with
    wide smearings -- the matrix is badly conditioned), WHITEN with C(t0)'s positive subspace: drop
    directions with eigenvalue < eps * max (noise), build C0^{-1/2} on the rest, and diagonalise the
    standard problem C0^{-1/2} C(t) C0^{-1/2}. The number of resolved states is the kept subspace
    dimension (< N if the operators are near-degenerate -- which is itself the answer: the basis
    only supports that many rungs)."""
    Nt, N, _ = Cmat.shape
    C0 = 0.5 * (Cmat[t0] + Cmat[t0].T)
    # the nucleon correlator can be uniformly NEGATIVE (a contraction-sign convention); then C0 is
    # negative-definite and sqrt(w0) is imaginary. Flip the whole matrix by its dominant-eigenvalue
    # sign so the signal subspace is positive, then whiten. The overall sign cancels in m_eff.
    wp = np.linalg.eigvalsh(C0)
    sgn = -1.0 if abs(wp.min()) > abs(wp.max()) else 1.0
    C0 = sgn * C0
    w0, V0 = np.linalg.eigh(C0)
    keep = w0 > eps * w0.max()
    lam = np.full((Nt, N), np.nan)
    if keep.sum() == 0:
        return lam
    Winv = V0[:, keep] / np.sqrt(w0[keep])                 # C0^{-1/2} on the kept subspace, (N, k)
    k = int(keep.sum())
    for t in range(Nt):
        Ct = sgn * 0.5 * (Cmat[t] + Cmat[t].T)
        A = Winv.T @ Ct @ Winv                            # k x k, well-conditioned
        ev = np.linalg.eigvalsh(A)
        lam[t, :k] = np.sort(ev)[::-1]
    return lam


def gevp_spectrum(raw, N, T, t0=2, tmin=None, tmax=None):
    """The excited torsiton spectrum from the variational (GEVP) correlator matrix.

    `raw` is the (n_rows, 5) array of rows [cfg, i, j, t, C_ij(t)] (measure_baryon_gevp output, tagged
    per config), i,j the N smearing operators. Solves the generalized eigenvalue problem
    C(t) v = lambda(t) C(t0) v on the config-averaged matrix: its N eigenvalues decay as
    lambda_n(t) ~ e^{-E_n (t-t0)}, so an N-operator basis resolves N states (ground .. highest).
    Returns the N masses (ascending) with jackknife-over-configs errors.

    The test of the three-generation thesis: how many clean nucleon levels appear, and whether a
    fourth fails to materialise. A state whose mass plateau is clean and stable is a real rung; one
    that is noise (huge error, no plateau) is the basis running out -- the tower's end."""
    raw = np.asarray(raw, dtype=float)
    cfg = raw[:, 0].astype(int)
    ii = raw[:, 1].astype(int)
    jj = raw[:, 2].astype(int)
    tt = raw[:, 3].astype(int)
    C = raw[:, 4]
    configs = np.unique(cfg)
    ncfg = len(configs)
    Nt = int(tt.max()) + 1

    def matrix(keep):
        M = np.zeros((Nt, N, N))
        cnt = np.zeros((Nt, N, N))
        np.add.at(M, (tt[keep], ii[keep], jj[keep]), C[keep])
        np.add.at(cnt, (tt[keep], ii[keep], jj[keep]), 1.0)
        return M / np.maximum(cnt, 1.0)

    lo = tmin if tmin is not None else max(t0 + 1, T // 6)
    hi = tmax if tmax is not None else max(lo + 1, T // 3)

    def fit(lam):
        with np.errstate(divide="ignore", invalid="ignore"):
            meff = np.log(lam[:-1] / lam[1:])                  # [Nt-1, N], E_n per state
        m = np.full(N, np.nan)
        for n in range(N):
            idx = np.arange(lo, min(hi, Nt - 1))
            idx = idx[np.isfinite(meff[idx, n])]
            if len(idx) >= 2:
                m[n] = float(np.mean(meff[idx, n]))
        return m, meff

    masses, meff = fit(_gevp_lambdas(matrix(np.ones(len(C), bool)), t0))
    masses = np.sort(masses)                                   # ascending: ground first
    jk = np.full((ncfg, N), np.nan)
    for a, cc in enumerate(configs):
        mj, _ = fit(_gevp_lambdas(matrix(cfg != cc), t0))
        jk[a] = np.sort(mj)
    err = np.full(N, np.nan)
    for n in range(N):
        v = jk[:, n][np.isfinite(jk[:, n])]
        if len(v) > 1:
            err[n] = float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))

    n_clean, count_verdict = _count_generations(masses, err)
    return {"masses": masses, "mass_err": err, "meff": meff, "t0": t0,
            "tmin": int(lo), "tmax": int(hi), "n_cfg": ncfg, "n_clean": n_clean,
            "verdict": count_verdict}


def _count_generations(masses, err, clean_frac=0.15):
    """Count the CLEAN radial rungs (= candidate generations) in a GEVP spectrum: a level is clean if
    it has a finite mass and a relative error below clean_frac, AND every lighter level is also clean
    (the tower must be contiguous from the ground state -- a clean level above a noisy gap is not a
    resolved rung). Returns (n_clean, verdict). The framework's live prediction is HOW MANY: three
    clean rungs and a noisy fourth postdicts the Standard Model's three generations; a clean fourth
    predicts a fourth (whose neutrino must be heavy to evade the LEP Z-width)."""
    masses = np.asarray(masses, dtype=float)
    err = np.asarray(err, dtype=float)
    n_clean = 0
    for n in range(len(masses)):
        ok = (np.isfinite(masses[n]) and np.isfinite(err[n]) and masses[n] > 0
              and err[n] < clean_frac * abs(masses[n]))
        if not ok:
            break
        n_clean += 1
    if n_clean <= 2:
        verdict = (f"{n_clean} clean rung(s) -- the basis/statistics resolve too few to test the "
                   f"generation count; add operators, configs, or a bigger box.")
    elif n_clean == 3:
        verdict = ("THREE clean rungs (and no clean fourth) -- consistent with exactly three "
                   "generations, a postdiction of the Standard Model's flavour count. Confirm the "
                   "fourth is genuinely noise (not just under-resolved) before claiming the cap.")
    else:
        verdict = (f"{n_clean} clean rungs -- a FOURTH (or higher) generation is resolved. If real, "
                   f"its neutrino must be heavy/sterile to evade the LEP Z-width (3 light neutrinos). "
                   f"Check it is a bound rung, not a two-particle scattering state.")
    return n_clean, verdict


# The lattice theory for the torsiton (L4). NOT the sextet -- that was a walking proxy for the
# electroweak S parameter (L3), a different question, and the wrong turn for the soliton (CLAUDE.md:
# "the sextet is DEAD"). The torsiton is the SU(3)-fundamental baryon (colour from the Pauli label),
# QCD-like; its mass and spectrum are standard baryon spectroscopy (measure_baryon).
CANDIDATE_THEORY = {
    "gauge_group": "SU(3)",
    "fermions": "fundamental Dirac quarks (colour forced by the Pauli label)",
    "object": "the torsiton = the SU(3)-fundamental baryon (three quarks, one of each colour)",
    "why": "the faithful realisation of the torsiton soliton; QCD-like, confining, chiral-breaking. "
           "The baryon ground-state mass (and excited spectrum) is standard lattice spectroscopy.",
    "gate_observable": "the nucleon (torsiton) effective-mass plateau m_N; the pion m_pi calibrates "
                       "the quark mass and probes the chiral condensate.",
    "first_step": "quenched valence on the L1 pure-gauge ensemble (measure_baryon), then dynamical "
                  "fundamental fermions near the chiral limit for the real spectrum.",
}
