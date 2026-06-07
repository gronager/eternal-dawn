"""The de-doubled (SLAC) chiral Dirac sea: free-field condensate structure vs Appendix C."""
import numpy as np

from cartasis_sims import dirac_sea as dsea


def test_slac_derivative_is_antisymmetric_and_chiral():
    D = dsea.slac_derivative(16, 0.3)
    assert np.allclose(D, -D.T)                     # anti-Hermitian: preserves chirality (no Wilson mass)
    assert np.allclose(np.diag(D), 0.0)
    # differentiates a smooth localised function in the interior (band-limited -> SLAC accurate)
    r = np.arange(16) * 0.3
    f = np.exp(-((r - 2.4) ** 2))
    dfdr = -2 * (r - 2.4) * f
    assert np.max(np.abs((D @ f - dfdr)[5:11])) < 5e-2


def test_free_sea_has_condensate_structure():
    # <qbar q> = -4 N_c M I0: negative, and |.| grows with M (Appendix C, eq C.3)
    vals = [dsea.free_sea_scalar_sum(M) for M in (0.5, 1.0, 1.5, 2.0)]
    assert all(v < 0 for v in vals)                 # condensate is negative
    assert vals[0] > vals[1] > vals[2] > vals[3]    # more negative with larger M
