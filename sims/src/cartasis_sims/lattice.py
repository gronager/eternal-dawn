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
    A = np.vstack([np.ones_like(x), x]).T
    AtW = A.T * w
    cov = np.linalg.inv(AtW @ A)
    c0, c1 = cov @ (AtW @ y)
    c0_err = float(np.sqrt(cov[0, 0]))
    rising = bool(c1 < 0)                               # s_T grows as m_pi2 falls toward chiral
    in_window = bool(lo <= c0 <= hi)
    lightest_s_T = float(y[np.argmin(x)])               # s_T at the lightest (most chiral) point
    any_in = bool(np.any((y >= lo) & (y <= hi)))        # is any MEASURED point already in the window?
    # a linear fit in m_pi^2 is a conservative LOWER bound when the rise accelerates toward chiral
    # (the lightest gaps widen); flag that so a c0 just-below-window is read honestly.
    accel = bool(len(x) >= 3 and rising and lightest_s_T > c0 + c0_err)

    span = float("nan")
    s_for_span = c0 if in_window else (lightest_s_T if any_in else c0)
    if in_window or any_in:
        from . import fermion_masses as fm
        lad = fm._ascending_ladder(float(np.clip(s_for_span, 0.30, 1.5)), well=fm.WELL)
        span = float(lad[-1] / lad[0])

    if in_window:
        verdict = (f"the bag GROWS into the productive window: chiral s_T = {c0:.2f}+/-{c0_err:.2f} r0 "
                   f"in [{lo},{hi}] -> lever span ~{span:.0f} vs observed 3477. CONSISTENT WITH "
                   f"DERIVED (valence chiral limit; unitary + 3-pt confirmation owed).")
    elif rising and c0 < lo and any_in:
        verdict = (f"ENCOURAGING: the linear chiral extrapolation sits just below the window "
                   f"(s_T={c0:.2f}+/-{c0_err:.2f} r0), but the lightest measured point is already IN it "
                   f"(s_T={lightest_s_T:.2f}, lever span ~{span:.0f})"
                   + (" and the rise ACCELERATES toward chiral (linear is a lower bound)" if accel else "")
                   + ". Firm it with more light-mass stats / the 3-pt condensate.")
    elif rising and c0 < lo:
        verdict = (f"the bag RISES toward chiral but extrapolates to s_T = {c0:.2f}+/-{c0_err:.2f} r0, "
                   f"still BELOW the window [{lo},{hi}] -- the one-body bag under-shoots; push to "
                   f"lighter mass and measure the 3-pt condensate <N|qbar q(r)|N>.")
    elif not rising:
        verdict = (f"the bag does NOT grow toward chiral (slope {c1:+.2f}); chiral s_T = {c0:.2f} r0. "
                   f"The one-body proxy under-shoots -- the hierarchy needs the 3-pt condensate.")
    else:
        verdict = (f"chiral s_T = {c0:.2f}+/-{c0_err:.2f} r0 ABOVE the window [{lo},{hi}] -- bag over-"
                   f"broad; the lever over-/under-shoots, recheck the scale and window.")

    return {"chiral_s_T": float(c0), "chiral_s_T_err": c0_err, "slope": float(c1), "rising": rising,
            "in_window": in_window, "span": span, "points": pts, "fit": (float(c0), float(c1)),
            "lightest_s_T": lightest_s_T, "any_in_window": any_in, "accelerating": accel,
            "verdict": verdict}


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
    w0, V0 = np.linalg.eigh(C0)
    keep = w0 > eps * w0.max()
    lam = np.full((Nt, N), np.nan)
    if keep.sum() == 0:
        return lam
    Winv = V0[:, keep] / np.sqrt(w0[keep])                 # C0^{-1/2} on the kept subspace, (N, k)
    k = int(keep.sum())
    for t in range(Nt):
        Ct = 0.5 * (Cmat[t] + Cmat[t].T)
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
    return {"masses": masses, "mass_err": err, "meff": meff, "t0": t0,
            "tmin": int(lo), "tmax": int(hi), "n_cfg": ncfg}


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
