r"""Deterministic well -> generation spectrum: measure the beast, don't fit it.

The Woods-Saxon bag (`dirac_woods_saxon`) welds the spectrum to a 2-parameter SHAPE. That shape
was always a stand-in for the one object the dynamics actually fixes: the dressed-quark self-energy,
i.e. the scalar mass function M(r) (the Fourier transform of the lattice Landau-gauge M(p^2)) and
its vector partner V(r). There is exactly ONE such pair for the theory at a given Lambda -- not a
family -- so once it is measured the spectrum has no shape freedom left.

This module is that pipeline. It takes ANY (r, M, V) -- a self-consistent mean-field well, a
parametric test well, or a tabulated LATTICE M(r) dropped in -- and returns, as OUTPUTS (never
inputs):

  * n_bound          -- how many radial levels the well actually binds (tests "exactly 3");
  * per level n:  E_n, node count, and two configurational-mass pieces with OPPOSITE rung-dependence
      m_local  = \int (G^2 - F^2) M(r) dr        (Eq. configmass: overlap with the local mass;
                                                  large for the spread-out rungs -> ground = lightest)
      m_core   = \int (G^2 + F^2)(M_vac - M) dr   (overlap with the condensate depletion / field-energy
                                                  source; large for the localized rung -> ground = heaviest)
    The two pieces are the quantitative form of the open "which rung is the electron?" fork.
  * the spans (heaviest/lightest) under each piece, vs the observed lepton/quark spans;
  * the sensitivity d ln(span)/d ln(shape) -- i.e. the lattice precision the hierarchy demands.

The relativistic first-order (G,F) Dirac is solved here by a SHOOTING integrator (RK4 on the well),
not a finite-difference matrix: the matrix solver (`dirac_soliton.dirac_levels`) suffers fermion
doublers -- grid-scale spurious modes that masquerade as radial excitations and so corrupt exactly
the generation levels we need. Shooting integrates the ODE, so the levels are clean and node-labelled
(the node theorem gives the radial quantum number n directly). This module adds that engine, the
shape-agnostic front end, the level bookkeeping, and the reporting.

Units are arbitrary (lengths in r0, masses in 1/r0): only ratios are physical, so c, hbar drop out.
"""
from __future__ import annotations

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz


# --- shooting engine: RK4 + inward/outward matching of the radial Dirac (no matrices, no doublers) -

def _deriv(kappa, E, ri, Mi, Vi, g, f):
    r"""RHS of  dG/dr = -(kappa/r)G + (E - V + M)F,   dF/dr = +(kappa/r)F - (E - V - M)G."""
    return (-(kappa / ri) * g + (E - Vi + Mi) * f,
            +(kappa / ri) * f - (E - Vi - Mi) * g)


def _rk4_step(kappa, E, r, M, V, i, j, g, f):
    """One RK4 step from grid index i to j (j=i+-1), interpolating M,V at the half point."""
    h = r[j] - r[i]
    rh = 0.5 * (r[i] + r[j]); Mh = 0.5 * (M[i] + M[j]); Vh = 0.5 * (V[i] + V[j])
    k1g, k1f = _deriv(kappa, E, r[i], M[i], V[i], g, f)
    k2g, k2f = _deriv(kappa, E, rh, Mh, Vh, g + 0.5 * h * k1g, f + 0.5 * h * k1f)
    k3g, k3f = _deriv(kappa, E, rh, Mh, Vh, g + 0.5 * h * k2g, f + 0.5 * h * k2f)
    k4g, k4f = _deriv(kappa, E, r[j], M[j], V[j], g + h * k3g, f + h * k3f)
    return (g + h / 6.0 * (k1g + 2 * k2g + 2 * k3g + k4g),
            f + h / 6.0 * (k1f + 2 * k2f + 2 * k3f + k4f))


def _solve_matched(E, r, M, V, kappa, im):
    r"""Integrate the regular solution outward (0 -> im) and the decaying solution inward
    (N-1 -> im), each L2-normalised over its own segment, and return (G, F, W). The decaying
    asymptotic start uses F/G = -sqrt((M_vac-E)/(M_vac+E)). The eigenvalue condition is the
    pole-free Wronskian  W(E) = G_out F_in - F_out G_in = 0 at the matchpoint (its sign change
    brackets a level; unlike F/G it has no poles where G_match -> 0). G,F is the spliced
    (G-continuous) wavefunction for node counting."""
    N = len(r)
    Go = np.empty(im + 1); Fo = np.empty(im + 1)
    Go[0] = r[0] ** abs(kappa); Fo[0] = 0.0
    for i in range(im):
        Go[i + 1], Fo[i + 1] = _rk4_step(kappa, E, r, M, V, i, i + 1, Go[i], Fo[i])
    Mvac = M[-1]
    ratio = -np.sqrt(max(Mvac - E, 0.0) / (Mvac + E))
    Gi = np.empty(N - im); Fi = np.empty(N - im)
    Gi[-1] = 1e-6; Fi[-1] = 1e-6 * ratio
    for k in range(N - 1, im, -1):
        kk = k - im
        Gi[kk - 1], Fi[kk - 1] = _rk4_step(kappa, E, r, M, V, k, k - 1, Gi[kk], Fi[kk])
    no = np.sqrt(max(_trapz(Go**2 + Fo**2, r[:im + 1]), 1e-300))
    ni = np.sqrt(max(_trapz(Gi**2 + Fi**2, r[im:]), 1e-300))
    Go /= no; Fo /= no; Gi /= ni; Fi /= ni
    W = Go[im] * Fi[0] - Fo[im] * Gi[0]                   # pole-free Wronskian
    scale = Go[im] / Gi[0] if Gi[0] != 0 else 0.0         # splice: G continuous at im
    G = np.concatenate([Go[:im], Gi * scale]); F = np.concatenate([Fo[:im], Fi * scale])
    return G, F, W


def _nodes(G):
    """Interior radial nodes of the large component = the radial quantum number n."""
    g = np.asarray(G)
    thr = 1e-6 * np.max(np.abs(g)) if g.size else 0.0
    sig = g[np.abs(g) > thr]
    if sig.size < 2:
        return 0
    return int(np.sum(np.diff(np.sign(sig)) != 0))


def _matchpoint(M):
    """Index of the well surface: where M(r) first rises through half its asymptotic value."""
    Mvac = M[-1]
    hit = np.where(M >= 0.5 * Mvac)[0]
    im = int(hit[0]) if hit.size else len(M) // 3
    return min(max(im, 3), len(M) - 3)


def dirac_shoot_levels(r, M, V, kappa=-1, nmax=12, nE=800):
    """Bound levels by inward/outward matching. Scan trial E in (0, M_vac), find sign changes of the
    matchpoint mismatch (each = one eigenvalue), bisect to pin E_n, normalise the spliced (G,F), and
    label by node count. Returns (list of (E,G,F) ascending in E, M_vac)."""
    M = np.asarray(M, float); V = np.asarray(V, float); r = np.asarray(r, float)
    Mvac = float(M[-1])
    im = _matchpoint(M)
    Es = np.linspace(1e-3 * Mvac, 0.9995 * Mvac, nE)
    D = np.array([_solve_matched(E, r, M, V, kappa, im)[2] for E in Es])
    levels = []
    for i in range(len(Es) - 1):
        d0, d1 = D[i], D[i + 1]
        if not (np.isfinite(d0) and np.isfinite(d1)) or d0 == 0 or d0 * d1 > 0:
            continue
        lo, hi = Es[i], Es[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            dm = _solve_matched(mid, r, M, V, kappa, im)[2]
            if d0 * dm <= 0:
                hi = mid
            else:
                lo = mid; d0 = dm
        E_n = 0.5 * (lo + hi)
        G, F, _ = _solve_matched(E_n, r, M, V, kappa, im)
        norm = np.sqrt(max(_trapz(G**2 + F**2, r), 1e-300))
        levels.append((E_n, G / norm, F / norm))
    levels.sort(key=lambda t: t[0])
    return levels[:nmax], Mvac

# observed spans for comparison (MeV)
LEPTONS = np.array([0.510999, 105.658, 1776.86])          # e, mu, tau
LEPTON_SPAN = float(LEPTONS[-1] / LEPTONS[0])             # ~3477
UP_QUARKS = np.array([2.16, 1270.0, 172_690.0])           # u, c, t (MS-bar-ish, illustrative)
DOWN_QUARKS = np.array([4.67, 93.4, 4180.0])              # d, s, b


def spectrum_from_well(r, M, V=None, kappa=-1, nmax=12):
    """Bound spectrum of the first-order scalar+vector Dirac in a GIVEN well.

    Parameters
    ----------
    r : (N,) radial grid (ascending, r[0] > 0).
    M : (N,) scalar mass function M(r) -- the dressed/constituent mass: ~M_vac outside, dipping in
        the chiral-restored core. (This IS the well; no shape is assumed.)
    V : (N,) or None -- vector self-energy V(r) (repulsive; enters as E -> E - V). None = scalar-only.
    kappa : Dirac angular quantum number (-1 = ground s_{1/2} channel).

    Returns
    -------
    dict with: 'n_bound', 'M_vac', and 'levels' = list (ascending E) of per-level dicts
    {'E','nodes','m_local','m_core','norm_in'}.
    """
    r = np.asarray(r, dtype=float)
    M = np.asarray(M, dtype=float)
    V = np.zeros_like(r) if V is None else np.asarray(V, dtype=float)
    lv, Mvac = dirac_shoot_levels(r, M, V, kappa=kappa, nmax=nmax)
    levels = []
    for E, G, F in lv:
        norm = _trapz(G**2 + F**2, r)
        m_local = _trapz((G**2 - F**2) * M, r)
        m_core = _trapz((G**2 + F**2) * (Mvac - M), r)
        levels.append({"E": float(E), "nodes": _nodes(G),
                       "m_local": float(m_local), "m_core": float(m_core),
                       "norm_in": float(norm)})
    levels.sort(key=lambda d: d["E"])
    return {"n_bound": len(levels), "M_vac": float(Mvac), "levels": levels}


def spans(res):
    """(span_local, span_core) = heaviest/lightest configurational mass under each piece, or None."""
    L = res["levels"]
    if len(L) < 2:
        return None, None
    ml = np.array([d["m_local"] for d in L])
    mc = np.array([d["m_core"] for d in L])
    sl = float(ml.max() / ml.min()) if ml.min() > 0 else float("inf")
    sc = float(mc.max() / mc.min()) if mc.min() > 0 else float("inf")
    return sl, sc


def well_from_mass_function(x, Mx, r, kind="r"):
    """Drop-in for a LATTICE-measured profile. Interpolate a tabulated mass function onto the
    solver grid r. kind='r': Mx is M(x) in position space (x=r). kind='p': Mx is M(p^2) and is
    transformed by a sine transform to M(r) first. Returns M(r) on the grid (constant-extrapolated
    outside the table). The vector V(r) is supplied the same way separately."""
    x = np.asarray(x, dtype=float)
    Mx = np.asarray(Mx, dtype=float)
    if kind == "p":
        # M(r) ~ (1/r) \int_0^\infty p M(p) sin(p r) dp / (2 pi^2) -- crude radial sine transform
        rr = np.asarray(r, dtype=float)
        Mr = np.array([np.trapz(x * Mx * np.sin(x * ri), x) / (2 * np.pi**2 * ri) for ri in rr])
        # normalise so the IR plateau matches M(p->0)
        Mr *= Mx[0] / max(Mr[np.argmin(np.abs(rr - rr.max()))], 1e-12)
        return Mr
    return np.interp(r, x, Mx, left=Mx[0], right=Mx[-1])


def shape_sensitivity(make_well, base, rel=0.02):
    """d ln(span)/d ln(s) about a 1-parameter family of wells. make_well(s)->(r,M,V); base=s0.
    Central difference on log span (m_local piece). Tells us the lattice precision the span needs."""
    s_lo, s_hi = base * (1 - rel), base * (1 + rel)
    sl_lo = spans(spectrum_from_well(*make_well(s_lo)))[0]
    sl_hi = spans(spectrum_from_well(*make_well(s_hi)))[0]
    if not sl_lo or not sl_hi or sl_lo <= 0 or sl_hi <= 0:
        return float("nan")
    return (np.log(sl_hi) - np.log(sl_lo)) / (np.log(s_hi) - np.log(s_lo))


def report(r, M, V=None, label="well"):
    """Print the honest read of a given well: level count, per-level E/nodes/masses, spans."""
    res = spectrum_from_well(r, M, V=V)
    print(f"[{label}]  M_vac = {res['M_vac']:.3f}   n_bound = {res['n_bound']}")
    if not res["levels"]:
        print("  (no bound levels)")
        return res
    print("   n  nodes   E/M_vac    m_local      m_core      norm_in")
    for k, d in enumerate(res["levels"]):
        print(f"  {k:>2}   {d['nodes']:>3}   {d['E']/res['M_vac']:>7.3f}   "
              f"{d['m_local']:>9.4f}   {d['m_core']:>9.4f}   {d['norm_in']:>6.3f}")
    sl, sc = spans(res)
    if sl is not None:
        print(f"  span(m_local, ground=lightest) = {sl:>8.1f}   "
              f"span(m_core, ground=heaviest) = {sc:>8.1f}")
        print(f"  observed: lepton span = {LEPTON_SPAN:.0f}  "
              f"(up-quark ~{UP_QUARKS[-1]/UP_QUARKS[0]:.0f}, down-quark ~{DOWN_QUARKS[-1]/DOWN_QUARKS[0]:.0f})")
    return res


# --- convenience wells for validation -------------------------------------------------

def woods_saxon_well(s_T, m_vac=0.9, a_frac=0.15, R_mult=12.0, N=900):
    """A WS test well in the (r, M, V) interface -- so the new pipeline can be cross-checked
    against `dirac_woods_saxon`. M(r)=m_vac/(1+exp(-(r-R0)/a)), R0=s_T (r0=1), V=0."""
    R0 = float(s_T)
    a = max(a_frac * R0, 1e-3)
    R = max(10.0, R_mult * R0)
    r = np.linspace(R / N, R, N)
    M = m_vac / (1.0 + np.exp(-(r - R0) / a))
    return r, M, np.zeros_like(r)


def well_from_bag_profile(prof, m_vac=0.9, r0_over_a=3.166, plateau=(4, 10), N=900, R_mult=4.0):
    """Turn a MEASURED lattice bag profile (measure_bag_profile PROF rows) into the scalar well M(r)
    -- the measured shape itself, no Woods-Saxon assumption. The gauge-invariant dressed-quark density
    rho(r) (peak=1 in the chiral-restored core, ->0 outside) maps to

        M(r) = m_vac * (1 - rho(r)/rho(0)),

    i.e. the mass melts where the quark sits and recovers to m_vac outside. Lengths are converted to
    r0 units via r0_over_a so the result is directly comparable to woods_saxon_well. `prof` is the
    (n,5) [cfg, r2, t, rho_sum, cnt] array (or any object lattice._bag_rho accepts). Returns (r,M,V=0).
    The vector well V(r) is a separate measurement (the vector self-energy / form factor)."""
    from . import lattice as lat
    r_lat, rho = lat._bag_rho(np.asarray(prof, dtype=float), plateau)
    r_phys = r_lat / r0_over_a                            # lattice spacings -> r0 units
    M_meas = m_vac * np.clip(1.0 - rho, 0.0, None)        # rho already normalised to peak=1
    R = max(R_mult * r_phys[-1], 10.0)
    r = np.linspace(R / N, R, N)
    M = np.interp(r, r_phys, M_meas, left=M_meas[0], right=m_vac)
    return r, M, np.zeros_like(r)


def report_bag_spectrum(prof, m_vac=0.9, r0_over_a=3.166, plateau=(4, 10), label="lattice bag"):
    """The (2) bridge end-to-end: a measured bag profile -> the scalar well -> the predicted
    generation spectrum (level count, node-labelled E_n, configurational-mass span). Run on the
    measure_bag_profile output to read the spectrum off the measured shape, no fitted well."""
    r, M, V = well_from_bag_profile(prof, m_vac=m_vac, r0_over_a=r0_over_a, plateau=plateau)
    return report(r, M, V=V, label=label)


def report_meanfield(g=4.0, v=1.0, lam=8.0, n_fermions=1, g_v=0.0, m_omega=2.0, N=240, R=10.0):
    """Feed the SELF-CONSISTENT mean-field well (dirac_soliton) into the pipeline and read its
    honest level count + span. This is the 'what does the un-fitted well actually predict?' test."""
    from . import dirac_soliton as ds
    out = ds.solve_dirac_soliton(g=g, v=v, lam=lam, n_fermions=n_fermions, g_v=g_v,
                                 m_omega=m_omega, N=N, R=R)
    r = out["r"]
    M = np.abs(g * out["sigma"])
    V = None
    return report(r, M, V=V, label=f"mean-field soliton (g={g}, g_v={g_v})")


def main():
    print("=== validation A: Woods-Saxon well through the shape-agnostic pipeline ===")
    for s in (2.0, 6.0, 8.0, 14.0):
        r, M, V = woods_saxon_well(s)
        report(r, M, V, label=f"WS s_T={s}")
        print()
    print("=== validation B: the un-fitted self-consistent mean-field well ===")
    report_meanfield(g=4.0, g_v=0.0)
    print()
    print("=== sensitivity: d ln(span)/d ln(s_T) for the WS family (lattice precision needed) ===")
    dlog = shape_sensitivity(lambda s: woods_saxon_well(s), base=8.0)
    print(f"  d ln(span_local)/d ln(s_T) ~ {dlog:.2f}  (|.|>>1 => the span is exponentially "
          f"sensitive to the measured shape)")


if __name__ == "__main__":
    main()
