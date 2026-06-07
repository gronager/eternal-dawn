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


def test_box_sea_reproduces_continuum_number():
    # with enough partial waves the box mode sum hits the continuum 3-momentum condensate (~3%,
    # the residual is finite-box). The number, not just the sign.
    import numpy as _np
    from cartasis_sims import dirac_sea as _ds
    M, Lam = 1.0, 4.0
    box = _ds.box_sea_condensate(M, Lam, kmax=30)
    cont = _ds.continuum_condensate(M, Lam)
    assert abs(box / cont - 1.0) < 0.06          # matches the continuum number


def test_soliton_sea_reinforces_bag_and_is_localized():
    # vacuum-subtracted sea condensate in a tanh bag: POSITIVE in the core (restores chiral symmetry
    # alongside the valence -> reinforces the bag) and LOCALISED (vanishes outside the bag)
    r = np.linspace(9.0 / 200, 9.0, 200)
    M = np.tanh(r / 1.0)                               # bag well, M_vac = 1
    sub = dsea.soliton_sea_condensate(M, 4.0, r, kmax=30)
    assert sub[np.argmin(np.abs(r - 0.7))] > 0.1      # positive (reinforcing) in the bag
    assert abs(sub[np.argmin(np.abs(r - 4.0))]) < 0.05  # localised: ~0 outside
