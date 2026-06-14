"""Shape-agnostic well -> spectrum pipeline: doubler-free, clean node labelling, deterministic."""
import numpy as np

from cartasis_sims import well_spectrum as ws


def test_clean_node_labelling_one_level_per_n():
    # a bag deep/wide enough for several radial levels must return EXACTLY one level per node count
    # n = 0, 1, 2, ... (the node theorem) -- the doubler-free property the matrix solver lacked.
    r, M, V = ws.woods_saxon_well(8.0, m_vac=0.9)
    res = ws.spectrum_from_well(r, M, V)
    nodes = [d["nodes"] for d in res["levels"]]
    assert res["n_bound"] >= 3
    assert nodes == list(range(len(nodes)))            # 0,1,2,... no gaps, no duplicates


def test_energies_ascending_and_in_gap():
    r, M, V = ws.woods_saxon_well(8.0, m_vac=0.9)
    res = ws.spectrum_from_well(r, M, V)
    E = [d["E"] for d in res["levels"]]
    assert E == sorted(E)
    assert all(0 < e < res["M_vac"] for e in E)        # genuine bound states in the mass gap


def test_level_count_grows_with_bag_size():
    # a bigger bag binds at least as many levels (monotone) -- the count is an honest output
    counts = [ws.spectrum_from_well(*ws.woods_saxon_well(s, m_vac=0.9))["n_bound"]
              for s in (4.0, 8.0, 14.0)]
    assert counts[0] <= counts[1] <= counts[2]
    assert counts[2] > counts[0]


def test_drop_in_tabulated_profile_matches_analytic():
    # the lattice drop-in (tabulate M(r), interpolate onto the solver grid) must reproduce the
    # spectrum of the analytic well it samples
    r0, M0, _ = ws.woods_saxon_well(8.0, m_vac=0.9)
    res0 = ws.spectrum_from_well(r0, M0)
    Mr = ws.well_from_mass_function(r0[::3], M0[::3], r0, kind="r")   # coarse table -> interp back
    res1 = ws.spectrum_from_well(r0, Mr)
    assert res1["n_bound"] == res0["n_bound"]
    e0 = [d["E"] for d in res0["levels"]]; e1 = [d["E"] for d in res1["levels"]]
    assert np.allclose(e0, e1, atol=0.02)


def _synth_bag(R0_lat=6.0, a=1.2, ncfg=2, tlo=4, thi=10, r2max=200):
    """A measure_bag_profile-format PROF array: dressed-quark density ~ Fermi bag of radius R0_lat."""
    rows = [[cfg, r2, t, 1.0 / (1.0 + np.exp((np.sqrt(r2) - R0_lat) / a)), 1.0]
            for cfg in range(1, ncfg + 1) for t in range(tlo, thi + 1) for r2 in range(0, r2max + 1)]
    return np.array(rows)


def test_lattice_bag_bridge_runs_and_binds():
    # the (2) bridge: a measured bag profile -> scalar well -> a clean bound spectrum
    prof = _synth_bag(R0_lat=6.0)
    r, M, V = ws.well_from_bag_profile(prof, m_vac=4.0)
    res = ws.spectrum_from_well(r, M, V)
    assert res["n_bound"] >= 1
    assert [d["nodes"] for d in res["levels"]] == list(range(res["n_bound"]))  # clean labelling
    assert M[0] < 0.5 * res["M_vac"]                   # melted core (mass off where the quark sits)
    assert abs(M[-1] - 4.0) < 1e-6                      # recovers to m_vac outside


def test_generation_count_set_by_bag_depth():
    # at a FIXED measured shape, the number of bound levels (= generations) grows with the depth
    # m_vac*r0: "exactly 3" is a constraint on the bag DEPTH, not only its sharpness
    prof = _synth_bag(R0_lat=6.0)
    counts = [ws.spectrum_from_well(*ws.well_from_bag_profile(prof, m_vac=mv))["n_bound"]
              for mv in (0.9, 4.0, 16.0)]
    assert counts[0] < counts[1] < counts[2]


def test_two_configmass_pieces_have_opposite_ordering():
    # the open "which rung is the electron" fork: m_local weighs spread-out rungs (ground lightest),
    # m_core weighs localized rungs (ground heaviest) -- they must order oppositely in n
    r, M, V = ws.woods_saxon_well(8.0, m_vac=0.9)
    res = ws.spectrum_from_well(r, M, V)
    ml = np.array([d["m_local"] for d in res["levels"]])
    mc = np.array([d["m_core"] for d in res["levels"]])
    assert ml[0] < ml[-1]                              # m_local: ground is lightest
    assert mc[0] > mc[-1]                              # m_core:  ground is heaviest
