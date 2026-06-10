"""Synthetic-data tests for the lattice analysis extractors (validated independent of any GPU run)."""

import numpy as np

from cartasis_sims import lattice as lat


def test_string_tension_recovered_from_synthetic_potential():
    # build a static potential with a KNOWN string tension and recover it
    r = np.linspace(1.0, 10.0, 20)
    sigma_true, alpha_true, c_true = 0.045, 0.30, 0.55
    V = c_true - alpha_true / r + sigma_true * r
    fit = lat.static_potential_cornell(r, V)
    assert np.isclose(fit["sigma"], sigma_true, rtol=1e-6)
    assert np.isclose(fit["alpha"], alpha_true, rtol=1e-6)
    assert fit["r0_sommer"] > 0          # confining -> a finite Sommer scale


def test_potential_from_wilson_loops_recovers_tension():
    # synthetic timelike loops W(R,T) = exp(-V(R) T) with a known Cornell V(R); recover sigma
    sigma_true, alpha_true, c_true = 0.05, 0.28, 0.6
    Rgrid = np.arange(1, 9)
    Tgrid = np.arange(1, 9)
    Rs, Ts, Ws = [], [], []
    for r in Rgrid:
        V = c_true - alpha_true / r + sigma_true * r
        for t in Tgrid:
            Rs.append(r); Ts.append(t); Ws.append(np.exp(-V * t))
    R, V = lat.effective_potential(np.array(Rs), np.array(Ts), np.array(Ws))
    fit = lat.static_potential_cornell(R, V)
    assert np.isclose(fit["sigma"], sigma_true, rtol=1e-6)
    assert np.isclose(fit["alpha"], alpha_true, rtol=1e-6)


def test_potential_linear_fit_recovers_tension():
    # same synthetic loops, extracted by the -ln W vs T slope fit (the robust extractor)
    sigma_true, alpha_true, c_true = 0.05, 0.28, 0.6
    Rs, Ts, Ws = [], [], []
    for r in np.arange(1, 9):
        V = c_true - alpha_true / r + sigma_true * r
        for t in np.arange(1, 9):
            Rs.append(r); Ts.append(t); Ws.append(np.exp(-V * t))
    R, V = lat.potential_from_loops_fit(np.array(Rs), np.array(Ts), np.array(Ws), tmin=1, tmax=6)
    fit = lat.static_potential_cornell(R, V)
    assert np.isclose(fit["sigma"], sigma_true, rtol=1e-6)


def test_effective_mass_table_plateaus_at_true_potential():
    # synthetic loops W(R,T)=exp(-V(R)T): V_eff(R,T)=ln[W(R,T)/W(R,T+1)] is a flat plateau at V(R)
    sigma_true, alpha_true, c_true = 0.05, 0.28, 0.6
    Rs, Ts, Ws = [], [], []
    Vtrue = {}
    for r in np.arange(1, 6):
        V = c_true - alpha_true / r + sigma_true * r
        Vtrue[int(r)] = V
        for t in np.arange(1, 9):
            Rs.append(r); Ts.append(t); Ws.append(np.exp(-V * t))
    table = lat.effective_mass_table(np.array(Rs), np.array(Ts), np.array(Ws))
    for r, (Tmid, veff) in table.items():
        assert np.allclose(veff, Vtrue[r], rtol=1e-6)   # flat plateau exactly at V(R)
        assert Tmid[0] == 1.0


def test_string_tension_jackknife_recovers_value_with_small_error():
    # many configs of clean synthetic loops (tiny per-config noise) -> sigma recovered, small error
    rng = np.random.default_rng(0)
    sigma_true, alpha_true, c_true = 0.157, 0.28, 0.6
    cfgs, Rs, Ts, Ws = [], [], [], []
    for ci in range(40):
        for r in np.arange(1, 7):
            V = c_true - alpha_true / r + sigma_true * r
            for t in np.arange(1, 7):
                w = np.exp(-V * t) * (1.0 + 0.01 * rng.standard_normal())
                cfgs.append(ci); Rs.append(r); Ts.append(t); Ws.append(w)
    out = lat.string_tension_jackknife(np.array(cfgs), np.array(Rs), np.array(Ts), np.array(Ws),
                                       tmin=2, tmax=5)
    assert out["n_cfg"] == 40
    assert np.isclose(out["sigma"], sigma_true, atol=0.02)
    assert out["sigma_err"] < 0.02                      # small error with clean, plentiful data
    assert abs(out["alpha"] - alpha_true) < 0.1


def test_rmax_cap_stabilises_jackknife_against_noisy_tail():
    # clean Cornell at R=1..6, then two GARBAGE large-R points (R=7,8) like the lattice L/2 tail.
    # capping the fit at rmax=6 must give a far smaller jackknife error than letting 7,8 in.
    rng = np.random.default_rng(1)
    sigma_true, alpha_true, c_true = 0.157, 0.28, 0.6
    cfgs, Rs, Ts, Ws = [], [], [], []
    for ci in range(30):
        for r in range(1, 9):
            V = c_true - alpha_true / r + sigma_true * r
            noise = 0.02
            if r >= 7:                                  # steep, KEPT-but-noisy tail (like V~L/2)
                V += 0.4 * (r - 6)                      # rises too fast -> stays monotonic, gets kept
                noise = 0.3                             # but with large per-config scatter
            for t in range(1, 7):
                w = np.exp(-V * t) * (1 + noise * rng.standard_normal())
                cfgs.append(ci); Rs.append(r); Ts.append(t); Ws.append(w)
    cfgs, Rs, Ts, Ws = map(np.array, (cfgs, Rs, Ts, Ws))
    capped = lat.string_tension_jackknife(cfgs, Rs, Ts, Ws, tmin=2, tmax=5, rmin=1, rmax=6)
    wild = lat.string_tension_jackknife(cfgs, Rs, Ts, Ws, tmin=2, tmax=5, rmin=1, rmax=8)
    assert max(capped["R_used"]) <= 6
    assert np.isclose(capped["sigma"], sigma_true, atol=0.03)     # capped fit recovers sigma
    assert capped["sigma_err"] < wild["sigma_err"]                # and is far more stable


def test_beta_calibration_interpolates_target_plaquette():
    # synthetic monotonic <P>(beta); recover the beta whose plaquette hits the target
    betas = np.array([5.6, 5.8, 6.0, 6.2])
    plaqs = np.array([0.524, 0.548, 0.567, 0.585])     # monotone rising
    out = lat.beta_from_plaquette(betas, plaqs, target=0.567)
    assert out["interpolated"]
    assert np.isclose(out["beta"], 6.0, atol=1e-6)
    # target above the scanned range -> extrapolate, flagged not-interpolated
    out2 = lat.beta_from_plaquette(betas, plaqs, target=0.60)
    assert not out2["interpolated"]
    assert out2["beta"] > 6.2


def test_screening_gives_zero_tension():
    # a pure Coulomb (screened) potential has sigma ~ 0 -> no area law
    r = np.linspace(1.0, 10.0, 20)
    V = 0.4 - 0.3 / r
    fit = lat.static_potential_cornell(r, V)
    assert abs(fit["sigma"]) < 1e-6


def test_deconfinement_peak_is_off_grid():
    # a susceptibility peaked between scan points -> beta_c interpolated, not snapped to grid
    beta = np.linspace(5.6, 6.0, 9)
    beta_c_true = 5.785
    chi = 1.0 / (1.0 + ((beta - beta_c_true) / 0.05) ** 2)
    out = lat.deconfinement_beta_c(beta, chi)
    assert np.isclose(out["beta_c"], beta_c_true, atol=0.01)
    assert beta_c_true not in beta        # genuinely between grid points


def test_anomalous_dimension_recovered_from_mode_number():
    # nu(M) = A * M^{4/(1+gamma)} with a known gamma -> recover gamma
    gamma_true = 0.35
    M = np.linspace(0.05, 0.5, 40)
    nu = 12.0 * M ** (4.0 / (1.0 + gamma_true))
    out = lat.anomalous_dimension_from_mode_number(M, nu)
    assert np.isclose(out["gamma_m"], gamma_true, rtol=1e-3)


def test_free_field_mode_number_gives_gamma_zero():
    # a free theory has nu(M) ~ M^4 (gamma = 0)
    M = np.linspace(0.05, 0.5, 40)
    nu = 3.0 * M ** 4
    out = lat.anomalous_dimension_from_mode_number(M, nu)
    assert abs(out["gamma_m"]) < 1e-6


def test_mode_number_counter_recovers_power_law():
    # eigenvalues drawn so the cumulative count nu(M) ~ M^4 -> the counter + fit recover gamma_m=0
    rng = np.random.default_rng(0)
    eig = rng.random(200000) ** 0.25            # CDF(M) ~ M^4 on [0,1]
    M = np.linspace(0.1, 0.6, 30)
    nu = lat.mode_number_from_eigenvalues(eig, M)
    assert np.all(np.diff(nu) >= 0)              # monotone cumulative count
    out = lat.anomalous_dimension_from_mode_number(M, nu, window=(0.15, 0.55))
    assert abs(out["gamma_m"]) < 0.1


def test_free_wilson_mode_number_recovers_gamma_zero():
    # the EXACT free Wilson-Dirac spectrum (no gauge field) must give gamma_m = 0, slope = 4.
    # needs a fine momentum grid (32^3+): on a coarse lattice the IR spectrum is gapped, not M^4
    # -- the same 'too coarse a grid' lesson as the string tension, now for the gate observable.
    M = np.linspace(0.10, 0.70, 60)
    nu = lat.free_wilson_mode_number(32, 32, M, m=0.0)
    out = lat.anomalous_dimension_from_mode_number(M, nu, window=(0.20, 0.50))
    assert abs(out["slope"] - 4.0) < 0.15       # nu ~ M^4
    assert abs(out["gamma_m"]) < 0.04           # gamma_m ~ 0


def test_free_wilson_too_coarse_has_no_scaling_window():
    # the negative control: on 16^3 the free spectrum is gapped, so the fit CANNOT recover gamma=0
    M = np.linspace(0.10, 0.70, 60)
    nu = lat.free_wilson_mode_number(16, 16, M, m=0.0)
    out = lat.anomalous_dimension_from_mode_number(M, nu, window=(0.20, 0.50))
    assert abs(out["gamma_m"]) > 0.1            # coarse grid -> no trustworthy gamma_m


def test_kpm_moments_recover_gamma_on_smooth_spectrum():
    # KPM (Chebyshev moments -> Jackson-damped step) must reproduce the mode number on a smooth
    # nu(M)~M^4 spectrum (gamma_m = 0). Validates moments -> nu(M) -> gamma_m end to end.
    rng = np.random.default_rng(0)
    lam = 8.0 * rng.random(300000) ** 0.25       # CDF(lambda) = (lambda/8)^4 -> nu ~ M^4
    M = np.linspace(1.0, 4.0, 30)
    win = (1.2, 3.6)
    mu = lat.chebyshev_moments_from_eigenvalues(lam, N=400, xmax=64.2)
    nu_kpm = lat.mode_number_from_chebyshev_moments(mu, M, xmax=64.2)
    nu_exact = lat.mode_number_from_eigenvalues(lam, M)
    g_kpm = lat.anomalous_dimension_from_mode_number(M, nu_kpm, window=win)["gamma_m"]
    g_exact = lat.anomalous_dimension_from_mode_number(M, nu_exact, window=win)["gamma_m"]
    assert abs(g_kpm - g_exact) < 0.02           # KPM tracks the exact mode number's slope
    assert abs(g_kpm) < 0.05                      # and recovers gamma_m ~ 0


def test_kpm_moments_match_exact_count_free_wilson():
    # on the exact free Wilson spectrum, moments -> nu(M) must reproduce the direct eigenvalue
    # count where the spectrum is DENSE (many modes). At small M the count is a coarse staircase
    # and KPM correctly returns its smooth trend -- a pointwise match there would be comparing a
    # smooth curve to a step function, so we validate in the dense regime instead.
    M = np.linspace(2.0, 4.5, 24)
    lam = lat._free_wilson_abs_lambda(24, 24)    # |lambda| list on a 24^4 free lattice
    mu = lat.chebyshev_moments_from_eigenvalues(lam, N=600, xmax=64.2)
    nu_kpm = lat.mode_number_from_chebyshev_moments(mu, M, xmax=64.2)
    nu_exact = lat.mode_number_from_eigenvalues(lam, M)
    rel = np.abs(nu_kpm - nu_exact) / nu_exact
    assert np.median(rel) < 0.01                 # KPM reproduces the exact count in the dense bulk
    assert rel.max() < 0.05


def test_gradient_flow_w0_recovered():
    # construct t^2<E>(t) so that W = t d/dt(t^2 E) crosses 0.3 at a known t -> recover w0
    t = np.linspace(0.01, 4.0, 400)
    # choose t^2 E = 0.15 * t  => W = t * 0.15 = 0.3 at t = 2 => w0 = sqrt(2)
    t2E = 0.15 * t
    w0 = lat.gradient_flow_w0(t, t2E, ref=0.3)
    assert np.isclose(w0, np.sqrt(2.0), rtol=1e-2)


def test_baryon_spectrum_recovers_synthetic_masses():
    # synthetic correlators C(t) = exp(-m t) with known masses, a few "configs" with small noise;
    # baryon_spectrum should recover m_pi and m_N from the effective-mass plateau
    import numpy as np
    from cartasis_sims import lattice as lat
    rng = np.random.default_rng(0)
    T = 32
    m_pi, m_N = 0.25, 0.70
    rows = []
    for cfg in range(1, 9):
        for t in range(T):
            cpi = np.exp(-m_pi * t) * (1 + 0.01 * rng.standard_normal())
            cn = np.exp(-m_N * t) * (1 + 0.01 * rng.standard_normal())
            rows.append([cfg, t, cpi, cn])
    res = lat.baryon_spectrum(np.array(rows), T=T)
    assert abs(res["pion"]["mass"] - m_pi) < 0.02
    assert abs(res["nucleon"]["mass"] - m_N) < 0.03
    assert res["pion"]["mass_err"] < 0.02 and res["nucleon"]["mass_err"] < 0.05
    assert res["nucleon"]["mass"] > res["pion"]["mass"]      # the baryon is heavier (bound, m_N>0)


def test_gevp_recovers_excited_states():
    # synthetic N-operator correlator matrix C_ij(t) = sum_n Z[n,i] Z[n,j] exp(-E_n t) with KNOWN
    # energies; the GEVP should recover the ground AND excited masses (the excited-torsiton test)
    import numpy as np
    from cartasis_sims import lattice as lat
    rng = np.random.default_rng(1)
    N, T = 3, 24
    E = np.array([0.45, 0.80, 1.20])                 # ground + two excited (the three rungs)
    Z = np.array([[1.0, 0.7, 0.4],                   # Z[state, operator]
                  [0.5, 1.1, 0.9],
                  [0.3, 0.6, 1.3]])
    rows = []
    for cfg in range(1, 13):
        for t in range(T):
            for i in range(N):
                for j in range(N):
                    c = sum(Z[n, i] * Z[n, j] * np.exp(-E[n] * t) for n in range(N))
                    c *= (1 + 0.003 * rng.standard_normal())
                    rows.append([cfg, i, j, t, c])
    res = lat.gevp_spectrum(np.array(rows), N=N, T=T, t0=2, tmin=3, tmax=8)
    assert abs(res["masses"][0] - 0.45) < 0.03        # ground
    assert abs(res["masses"][1] - 0.80) < 0.06        # first excited (2nd generation)
    assert abs(res["masses"][2] - 1.20) < 0.15        # second excited (3rd generation)
    assert res["masses"][0] < res["masses"][1] < res["masses"][2]


def test_nodal_gevp_recovers_three_generations():
    # The NODAL-basis test: the covariant-Laplacian polynomial sources (0,1,2 radial nodes) give the
    # GEVP an N=3 cross-correlator matrix C_ij(t) = sum_n Z[n,i] Z[n,j] exp(-E_n t). With three
    # well-separated, cleanly-overlapping levels, gevp_spectrum must recover all three eigen-energies
    # AND _count_generations must report exactly 3 clean rungs (the three-generation postdiction).
    # This is the Python-side gate on the nodal pipeline; the C++ nodal smearing only changes the
    # OVERLAPS Z (giving real 1-/2-node sensitivity), not the on-disk [cfg,i,j,t,C] format this reads.
    import numpy as np
    from cartasis_sims import lattice as lat
    rng = np.random.default_rng(7)
    N, T = 3, 24
    E = np.array([0.40, 0.75, 1.15])                 # three KNOWN eigen-energies: gen I, II, III
    # node-bearing overlaps: operator k (k nodes) couples most strongly to state k, with controlled
    # off-diagonal mixing (a well-conditioned, non-collinear basis -- what the nodes buy us).
    Z = np.array([[1.0, 0.45, 0.20],                 # Z[state, operator]
                  [0.50, 1.0, 0.55],
                  [0.25, 0.50, 1.0]])
    rows = []
    for cfg in range(1, 13):
        for t in range(T):
            for i in range(N):
                for j in range(N):
                    c = sum(Z[n, i] * Z[n, j] * np.exp(-E[n] * t) for n in range(N))
                    c *= (1 + 0.002 * rng.standard_normal())
                    rows.append([cfg, i, j, t, c])
    res = lat.gevp_spectrum(np.array(rows), N=N, T=T, t0=2, tmin=3, tmax=8)
    assert abs(res["masses"][0] - E[0]) < 0.03        # ground (gen I)
    assert abs(res["masses"][1] - E[1]) < 0.06        # first excited (gen II)
    assert abs(res["masses"][2] - E[2]) < 0.15        # second excited (gen III)
    assert res["masses"][0] < res["masses"][1] < res["masses"][2]
    # the recovered tower is three clean, contiguous rungs -> three generations
    n_clean, verdict = lat._count_generations(res["masses"], res["mass_err"])
    assert n_clean == 3
    assert res["n_clean"] == 3


def _synthetic_bag(R0_true, a_true=0.5, L=12, peak=1.0, ncfg=5, plateau=(4, 10), seed=0):
    """measure_bag_profile-style PROF rows [cfg, r2, t, rho_sum, count] for a known Fermi bag
    rho_norm(r)=1/(1+exp((r-R0)/a)), with the exact lattice r^2 shell counts (nearest image)."""
    from collections import defaultdict
    cnt = defaultdict(int)
    for x in range(L):
        dx = min(x, L - x)
        for y in range(L):
            dy = min(y, L - y)
            for z in range(L):
                dz = min(z, L - z)
                cnt[dx * dx + dy * dy + dz * dz] += 1
    fermi = lambda r: 1.0 / (1.0 + np.exp((r - R0_true) / a_true))
    rng = np.random.default_rng(seed)
    rows = []
    for c in range(1, ncfg + 1):
        for t in range(plateau[0], plateau[1] + 1):
            for r2, k in cnt.items():
                rho_site = peak * fermi(np.sqrt(r2)) * (1.0 + 0.02 * rng.standard_normal())
                rows.append([c, r2, t, rho_site * k, k])     # rho_sum = per-site density x count
    return np.array(rows, dtype=float)


def test_bag_profile_recovers_half_density_radius():
    # a known bag of half-density radius R0=2.5 (lattice units) is recovered from rho(r)
    prof = _synthetic_bag(R0_true=2.5)
    res = lat.bag_profile(prof, T=24, plateau=(4, 10), r0_over_a=3.166)
    assert abs(res["R0"] - 2.5) < 0.3
    assert np.isclose(res["s_T"], res["R0"] / 3.166)
    assert res["span"] > 0 and res["n_cfg"] == 5


def test_bag_profile_sharper_bag_gives_larger_span():
    # the lever direction: a SMALLER (sharper) bag -> larger generation span
    sharp = lat.bag_profile(_synthetic_bag(R0_true=1.5), T=24, r0_over_a=3.166)
    broad = lat.bag_profile(_synthetic_bag(R0_true=3.5), T=24, r0_over_a=3.166)
    assert sharp["s_T"] < broad["s_T"]
    assert sharp["span"] > broad["span"]
    assert "SHARP" in sharp["verdict"] or sharp["span"] > broad["span"]


def test_correlator_mass_recovers_known_decay():
    # C(t) = A exp(-m t) across configs; recover m from the effective-mass plateau
    m_true, A, T, ncfg = 0.42, 3.0, 24, 6
    rng = np.random.default_rng(1)
    rows = []
    for c in range(1, ncfg + 1):
        for t in range(T):
            rows.append([c, t, A * np.exp(-m_true * t) * (1 + 0.01 * rng.standard_normal())])
    res = lat.correlator_mass(np.array(rows), T=T, tmin=4, tmax=10)
    assert abs(res["mass"] - m_true) < 0.02
    assert res["n_cfg"] == ncfg


def test_bag_chiral_trend_rises_into_window():
    # s_T grows as m_pi^2 -> 0 and extrapolates INTO the productive window [0.43, 0.70]
    mpi2 = np.array([1.7, 1.0, 0.5, 0.2])
    s_T = 0.58 - 0.18 * mpi2                     # intercept 0.58 (in window), negative slope = rising
    pts = list(zip(mpi2, s_T, 0.01 * np.ones_like(s_T)))
    tr = lat.bag_chiral_trend(pts)
    assert tr["rising"] is True
    assert 0.43 <= tr["chiral_s_T"] <= 0.70 and tr["in_window"]
    assert tr["span"] > 0 and "DERIVED" in tr["verdict"]


def test_bag_chiral_trend_undershoots_below_window():
    # rises toward chiral but extrapolates BELOW the window -> honest "under-shoots"
    mpi2 = np.array([1.7, 1.0, 0.5, 0.2])
    s_T = 0.38 - 0.07 * mpi2                     # intercept 0.38 < 0.43
    pts = list(zip(mpi2, s_T, 0.01 * np.ones_like(s_T)))
    tr = lat.bag_chiral_trend(pts)
    assert tr["rising"] is True and not tr["in_window"]
    assert tr["chiral_s_T"] < 0.43 and "under-shoots" in tr["verdict"]


def test_bag_chiral_trend_single_point_no_extrapolation():
    tr = lat.bag_chiral_trend([(1.7, 0.27, 0.01)])
    assert "single mass" in tr["verdict"] and np.isnan(tr["chiral_s_T_err"])


def _synthetic_3pt(R0_true, T=24, t_snk=12, cn=5.0, gS=0.8, ncfg=6, selfcheck="pass", seed=2):
    """measure_condensate_3pt-style rows: a known Fermi condensate bag G3(r) over a tau plateau, a
    2pt C_N(t_snk)=cn, a flat sum rule SR=gS*cn, and a CHK pair (recon vs C_N)."""
    from collections import defaultdict
    L = 12
    cnt = defaultdict(int)
    for x in range(L):
        dx = min(x, L - x)
        for y in range(L):
            dy = min(y, L - y)
            for z in range(L):
                dz = min(z, L - z)
                cnt[dx * dx + dy * dy + dz * dz] += 1
    fermi = lambda r: 1.0 / (1.0 + np.exp((r - R0_true) / 0.5))
    rng = np.random.default_rng(seed)
    chk, c2, sr, p3 = [], [], [], []
    recon = cn if selfcheck == "pass" else 0.3 * cn       # mismatch => self-check fails
    for c in range(1, ncfg + 1):
        chk.append([c, recon * (1 + 0.001 * rng.standard_normal()), cn])
        for t in range(T):
            c2.append([c, t, cn * np.exp(-0.0 * t)])      # flat enough near the sink for the test
            sr.append([c, t, gS * cn * (1 + 0.01 * rng.standard_normal())])
            for r2, k in cnt.items():
                p3.append([c, r2, t, gS * fermi(np.sqrt(r2)) * k, k])
    return (np.array(p3), np.array(c2), np.array(sr), np.array(chk))


def test_condensate_3pt_selfcheck_passes_and_recovers_bag():
    p3, c2, sr, chk = _synthetic_3pt(R0_true=1.4)
    res = lat.condensate_3pt(p3, c2, sr, chk, T=24, t_snk=12, tau_window=(4, 8), r0_over_a=3.166)
    assert res["self_check_ok"] is True
    assert abs(res["g_S"] - 0.8) < 0.05
    assert abs(res["R0"] - 1.4) < 0.3 and np.isclose(res["s_T"], res["R0"] / 3.166)
    assert "PASSED" in res["verdict"]


def test_condensate_3pt_selfcheck_failure_is_flagged():
    p3, c2, sr, chk = _synthetic_3pt(R0_true=1.4, selfcheck="fail")
    res = lat.condensate_3pt(p3, c2, sr, chk, T=24, t_snk=12, tau_window=(4, 8))
    assert res["self_check_ok"] is False
    assert "SELF-CHECK FAILED" in res["verdict"]


def test_condensate_3pt_far_sink_is_flagged():
    # self-check passes (contraction correct) but the sink sits past the 2pt node (C_N<0) -> flagged
    p3, c2, sr, chk = _synthetic_3pt(R0_true=1.4)
    c2 = c2.copy()
    c2[c2[:, 1] >= 10, 2] *= -1.0                 # nucleon 2pt flips sign at t=10 (backward state)
    chk = chk.copy()
    chk[:, 1] = -1.0; chk[:, 2] = -1.0            # recon==C_N(t_snk)<0: self-check ok, sink bad
    res = lat.condensate_3pt(p3, c2, sr, chk, T=24, t_snk=12, tau_window=(4, 8))
    assert res["self_check_ok"] is True and res["sink_ok"] is False
    assert res["node_t"] == 10 and "PAST THE PLATEAU" in res["verdict"]


def test_condensate_3pt_negative_convention_sink_ok():
    # a globally NEGATIVE nucleon 2pt (sign convention) with a clean plateau must NOT be flagged
    p3, c2, sr, chk = _synthetic_3pt(R0_true=1.4)
    c2 = c2.copy(); c2[:, 2] *= -1.0                  # flip the whole correlator negative
    chk = chk.copy(); chk[:, 1] *= -1.0; chk[:, 2] *= -1.0
    res = lat.condensate_3pt(p3, c2, sr, chk, T=24, t_snk=8, tau_window=(3, 6))
    assert res["self_check_ok"] is True and res["sink_ok"] is True   # sign-robust: not a node


def _synthetic_gevp(masses, T=24, N=3, t0=2, amp=None, sign=1.0, seed=3):
    """An N-operator GEVP correlator matrix C_ij(t) = sum_n A_in A_jn e^{-m_n t} for known masses,
    scaled by `sign` (to test the globally-negative convention). Rows [cfg, i, j, t, C]."""
    rng = np.random.default_rng(seed)
    masses = np.asarray(masses, float)
    A = np.eye(N) + 0.3 * rng.standard_normal((N, N)) if amp is None else amp
    rows = []
    for cfg in range(1, 9):
        for t in range(T):
            C = sign * (A * np.exp(-masses * t)[None, :]) @ A.T
            for i in range(N):
                for j in range(N):
                    rows.append([cfg, i, j, t, C[i, j] * (1 + 0.005 * rng.standard_normal())])
    return np.array(rows)


def test_gevp_sign_robust_recovers_masses_when_negative():
    # a globally NEGATIVE correlator matrix (sign convention) must give the SAME answer as positive
    true = [0.5, 0.9, 1.4]
    res_pos = lat.gevp_spectrum(_synthetic_gevp(true, sign=+1.0), N=3, T=24, t0=2, tmin=3, tmax=8)
    res_neg = lat.gevp_spectrum(_synthetic_gevp(true, sign=-1.0), N=3, T=24, t0=2, tmin=3, tmax=8)
    assert np.allclose(res_pos["masses"], res_neg["masses"], equal_nan=True)   # the point: sign-robust
    for n in range(2):                                   # the two cleanly-resolved levels
        assert abs(res_pos["masses"][n] - true[n]) < 0.05


def test_generation_count_three_vs_four():
    # three clean rungs -> "three generations"; a clean fourth -> "FOURTH ... generation"
    n3, v3 = lat._count_generations([0.5, 0.9, 1.4, 2.0], [0.01, 0.02, 0.05, 0.9])
    assert n3 == 3 and "THREE clean" in v3
    n4, v4 = lat._count_generations([0.5, 0.9, 1.4, 2.0], [0.01, 0.02, 0.05, 0.1])
    assert n4 == 4 and "FOURTH" in v4
    # a noisy second breaks the contiguous tower at one
    n1, _ = lat._count_generations([0.5, 0.9, 1.4], [0.01, 0.9, 0.05])
    assert n1 == 1


def _synthetic_3pt_twostate(R0_gs=1.4, R0_ex=0.7, dE=0.6, amp_ex=0.8, t_snk=12, T=24, L=12, ncfg=6, seed=5):
    """3-pt rows where G3(r,tau) = M_gs(r) + amp_ex*M_ex(r)*[e^{-dE tau}+e^{-dE(t_snk-tau)}], a known
    ground-state Fermi bag (R0_gs) plus a more-compact excited state (R0_ex<R0_gs). The tau-plateau
    average is contaminated; the two-state fit must recover R0_gs."""
    from collections import defaultdict
    cnt = defaultdict(int)
    for x in range(L):
        dx = min(x, L - x)
        for y in range(L):
            dy = min(y, L - y)
            for z in range(L):
                dz = min(z, L - z)
                cnt[dx * dx + dy * dy + dz * dz] += 1
    fg = lambda r, R: 1.0 / (1.0 + np.exp((r - R) / 0.5))
    g = lambda tau: np.exp(-dE * tau) + np.exp(-dE * (t_snk - tau))
    rng = np.random.default_rng(seed)
    cn = -1.0                                            # negative convention (sign-robust path)
    p3, sr, c2, chk = [], [], [], []
    for c in range(1, ncfg + 1):
        chk.append([c, cn, cn])
        for t in range(T):
            SR = 0.0
            for r2, k in cnt.items():
                val = fg(np.sqrt(r2), R0_gs) + amp_ex * fg(np.sqrt(r2), R0_ex) * g(t)
                SR += val * k
                p3.append([c, r2, t, val * k * (1 + 0.004 * rng.standard_normal()), k])
            sr.append([c, t, SR]); c2.append([c, t, cn])
    return np.array(p3), np.array(c2), np.array(sr), np.array(chk)


def test_twostate_fit_recovers_groundstate_shape():
    p3, c2, sr, chk = _synthetic_3pt_twostate(R0_gs=1.4, R0_ex=0.7)
    gs = lat._twostate_profile(p3, sr, t_snk=12, r0_over_a=3.166)
    assert gs is not None and abs(gs["R0"] - 1.4) < 0.3      # ground state recovered, not the mix
    assert abs(gs["dE"] - 0.6) < 0.2                         # the gap comes out right too


def test_condensate_3pt_reports_groundstate_bag():
    p3, c2, sr, chk = _synthetic_3pt_twostate(R0_gs=1.4, R0_ex=0.7)
    res = lat.condensate_3pt(p3, c2, sr, chk, T=24, t_snk=12, tau_window=(4, 8), r0_over_a=3.166)
    assert res["self_check_ok"] and abs(res["gs_R0"] - 1.4) < 0.3
    assert "GROUND-STATE" in res["verdict"]


def _synthetic_pp(m_pi, G_PP, T, ncfg=8, noise=0.005, seed=7):
    """measure_decay_constant-style PP rows [cfg, t, C_PP(t)] for a known two-sided cosh
    C_PP(t) = (G_PP/(2 m_pi)) (e^{-m_pi t} + e^{-m_pi (T-t)}), with small per-config noise."""
    rng = np.random.default_rng(seed)
    tt = np.arange(T)
    shape = (G_PP / (2.0 * m_pi)) * (np.exp(-m_pi * tt) + np.exp(-m_pi * (T - tt)))
    rows = []
    for c in range(1, ncfg + 1):
        for t in range(T):
            rows.append([c, t, shape[t] * (1.0 + noise * rng.standard_normal())])
    return np.array(rows, dtype=float)


def test_decay_constant_recovers_mpi_and_fpi():
    # synthetic cosh PP correlator with KNOWN m_pi and G_PP; decay_constant must recover both,
    # and f_pi*a must equal the PCAC formula (2 m_q / m_pi^2) sqrt(G_PP).
    m_pi_true, G_PP_true, m_q, T = 0.30, 0.8, 0.05, 32
    pp = _synthetic_pp(m_pi_true, G_PP_true, T)
    res = lat.decay_constant(pp, m_q=m_q, T=T, tmin=4, tmax=14)
    assert abs(res["m_pi"] - m_pi_true) < 0.02
    assert abs(res["G_PP"] - G_PP_true) < 0.05
    f_expected = (2.0 * m_q / m_pi_true**2) * np.sqrt(G_PP_true)
    assert abs(res["f_pi"] - f_expected) < 0.05 * f_expected
    assert np.isclose(res["v"], res["f_pi"])           # the framework's v = f_pi identification
    assert res["m_pi_err"] < 0.01 and res["f_pi_err"] < 0.05 * f_expected
    assert res["n_cfg"] == 8


def test_sigma_2piv2_check_flags_satisfied_and_violated():
    # construct (sigma, f_pi) pairs that exactly satisfy / grossly violate sigma = 2 pi v^2.
    f_pi_a = 0.10
    sqrt_sigma_sat = np.sqrt(2.0 * np.pi * f_pi_a**2)   # sigma = 2 pi f_pi^2 exactly -> ratio 1
    sat = lat.sigma_2piv2_check(sqrt_sigma_sat, f_pi_a)
    assert np.isclose(sat["ratio"], 1.0, atol=1e-9)
    assert sat["consistent"] and "CONSISTENT" in sat["verdict"]
    # gross violation: sigma 10x too large -> ratio ~ 10, outside the band
    viol = lat.sigma_2piv2_check(sqrt_sigma_sat * np.sqrt(10.0), f_pi_a)
    assert viol["ratio"] > 5.0 and not viol["consistent"]
    assert "STRAINED" in viol["verdict"]


def test_sigma_2piv2_check_within_errors():
    # a central value off the band but bracketing 1 within the propagated error -> within_err True
    f_pi_a = 0.10
    sqrt_sigma_sat = np.sqrt(2.0 * np.pi * f_pi_a**2)
    res = lat.sigma_2piv2_check(sqrt_sigma_sat * 1.3, f_pi_a,
                                sqrt_sigma_a_err=0.4 * sqrt_sigma_sat, f_pi_a_err=0.0)
    assert res["within_err"] and np.isfinite(res["ratio_err"])
