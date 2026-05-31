"""Tests for the extruder (entropy budget of the bounce)."""

import math

import numpy as np

from cartasis_sims import extruder as ex
from cartasis_sims import cve_filter as cve


def test_entropy_integral_is_order_unity():
    # J = int H^2 a dtau is a finite O(1) pure number for the radiation bounce
    res = ex.analyze(w=1.0 / 3.0)
    assert 0.5 < res.J < 5.0


def test_adiabaticity_equals_cve_small_parameter():
    # the prefactor Omega/T is exactly the CVE's omega/T (maximal-spin) parameter
    res = ex.analyze(w=1.0 / 3.0)
    assert math.isclose(res.adiabaticity, cve.vorticity_over_T(spin_fraction=1.0),
                        rel_tol=1e-9)


def test_conformal_fluid_conserves_entropy():
    # zeta_tilde = 0 (conformal radiation) -> no entropy generated -> D = 1
    assert ex.dilution_factor(0.0) == 1.0
    assert math.isclose(ex.inherited_eta(0.0), cve.ETA_OBS, rel_tol=1e-12)


def test_physical_viscosity_is_adiabatic_to_1e10():
    # even O(1) bulk viscosity leaves D-1 below ~1e-9: inheritance essentially exact
    res = ex.analyze(w=1.0 / 3.0)
    assert ex.dilution_factor(1.0, res) - 1.0 < 1e-9
    assert math.isclose(ex.inherited_eta(1.0, res=res), cve.ETA_OBS, rel_tol=1e-8)


def test_horizon_saturation_needs_absurd_viscosity():
    # reaching D ~ S_BH/S_rad would require zeta_tilde ~ 1e11-1e12, unphysical
    res = ex.analyze(w=1.0 / 3.0)
    D_hor = cve.dilution_horizon(1.0e53)
    zt_needed = math.log(D_hor) / (9.0 * res.adiabaticity * res.J)
    assert zt_needed > 1e10


def test_ln_dilution_linear_in_viscosity():
    res = ex.analyze(w=1.0 / 3.0)
    assert math.isclose(ex.ln_dilution(2.0, res), 2.0 * ex.ln_dilution(1.0, res),
                        rel_tol=1e-9)
