"""Tests for gravity-scaled coordinates (the conformal map)."""

import math

from cartasis_sims import gravity_scale as gs


def test_reference_density_cancels_in_ratios():
    # R_gs ratio depends only on the mass ratio, not on rho_ref
    r1 = gs.gravity_scaled_radius(1e50, rho_ref=1000.0)
    r2 = gs.gravity_scaled_radius(1e40, rho_ref=1000.0)
    r1b = gs.gravity_scaled_radius(1e50, rho_ref=1.0)
    r2b = gs.gravity_scaled_radius(1e40, rho_ref=1.0)
    assert math.isclose(r1 / r2, r1b / r2b, rel_tol=1e-9)
    assert math.isclose(r1 / r2, (1e50 / 1e40) ** (1 / 3), rel_tol=1e-9)


def test_ogu_is_not_magnified_vs_void():
    # a horizon-mass hole has ~ the cosmic/void density -> magnification ~ 1
    mu = gs.ogu_void_magnification(9.2e52)
    assert 0.8 < mu < 1.5                      # ~1.1: gravity-scaled = true at OGU level


def test_horizon_mass_hole_has_cosmic_density():
    # the R_s = R_H identity: an OGU's mean density ~ the void density
    rho = gs.black_hole_mean_density(9.2e52)
    assert math.isclose(rho, gs.RHO_LAMBDA, rel_tol=0.5)


def test_astrophysical_holes_are_hugely_magnified():
    # stellar / supermassive holes are far denser than their universe -> big mu
    assert gs.bh_universe_magnification(2e31) > 1e12     # stellar
    assert gs.bh_universe_magnification(2e39) > 1e8      # SMBH
    # and denser holes magnify more (mean density ~ 1/M^2)
    assert gs.bh_universe_magnification(2e31) > gs.bh_universe_magnification(2e39)


def test_bh_mean_density_falls_as_inverse_mass_squared():
    assert math.isclose(
        gs.black_hole_mean_density(1e40) / gs.black_hole_mean_density(1e41),
        100.0, rel_tol=1e-6)


def test_ogu_separation_over_diameter_is_absurdly_large():
    # the one emptiness number: sep/diam ~ 10^(2.7e83) for I ~ 2.5e84
    from cartasis_sims import birth_rate as br
    log10_ratio = gs.ogu_separation_over_diameter_log10(br.seed_instanton_action())
    assert 1e83 < log10_ratio < 1e84          # exponent ~2.7e83 (dwarfs 1e80 atoms)
