"""Tests for the Oppenheimer-Snyder torsion-bounce collapse."""

import numpy as np

from cartasis_sims import os_collapse as osc
from cartasis_sims import constants as k


def test_bounce_is_deep_inside_horizon():
    # For the observable-universe mass the bounce areal radius is metres while
    # the horizon is ~10^26 m: R_min/r_s ~ 1e-26.
    ratio = osc.rmin_over_rs(9.246e52, rho_C=1.0e50)
    assert ratio < 1e-20
    assert osc.r_min(9.246e52, 1.0e50) < 100.0       # a few metres
    assert osc.r_schwarzschild(9.246e52) > 1e25      # ~Hubble radius


def test_bounce_volume_matches_draft_scale():
    # ~10^3 m^3 at the bounce (the draft's "half an Olympic pool").
    import math
    R = osc.r_min(9.246e52, 1.0e50)
    vol = (4.0 / 3.0) * math.pi * R**3
    assert 1e2 < vol < 1e4


def test_exterior_time_diverges_at_horizon():
    R, t, z = osc.exterior_geodesic(R0_over_rs=3.0)
    # Schwarzschild time grows without bound as R -> r_s (frozen star)
    assert t[-1] > 5.0
    assert t[-1] == t.max()
    # redshift blows up at the horizon
    assert z[-1] > 5.0


def test_interior_bounces_not_singular():
    s = osc.simulate_os(a_init=8.0, rs_units=3.0)
    assert s.R_int.min() > 0.9          # bounces at R_min = 1 (units)
    assert s.R_os.min() < 0.1           # GR dives toward 0
    # the bounce sits below the (illustrative) horizon
    assert s.R_int.min() < s.r_s_units
