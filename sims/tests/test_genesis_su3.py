"""Horizon 3b: the Skyrme baryon octet+decuplet from one soliton (N, Delta, hyperons)."""
import numpy as np
from cartasis_sims import genesis_su3 as su3


def test_hedgehog_profile_solves():
    x, F, M_red = su3.hedgehog_profile()
    assert abs(F[0] - np.pi) < 1e-6 and abs(F[-1]) < 1e-3   # F(0)=pi -> F(inf)=0
    assert M_red > 0 and su3.reduced_inertia(x, F) > 0       # finite soliton mass and inertia


def test_baryon_spectrum_predicts_hyperons():
    res = su3.baryon_spectrum()
    m, obs = res["mass"], su3.OBS
    # N, Delta, Omega are the calibration -> exact; the rest are predictions
    for cal in ("N", "Delta", "Omega"):
        assert abs(m[cal] - obs[cal]) < 1.0
    # the predicted hyperons land within Skyrme-model accuracy
    for pred, tol in [("Sigma*", 0.04), ("Xi*", 0.04), ("Lambda", 0.06),
                      ("Sigma", 0.10), ("Xi", 0.08)]:
        assert abs(m[pred] - obs[pred]) / obs[pred] < tol
    # the decuplet EQUAL SPACING (Delta, Sigma*, Xi*, Omega) -- a clean prediction
    dec = [m["Delta"], m["Sigma*"], m["Xi*"], m["Omega"]]
    gaps = np.diff(dec)
    assert np.allclose(gaps, gaps.mean(), rtol=1e-6)        # exactly equal by construction (the model)
    assert 120 < gaps.mean() < 170                          # ~147 MeV, the observed spacing


def test_calibration_gives_sensible_scales():
    res = su3.baryon_spectrum()
    s = res["scales"]
    assert 60 < s["f_pi"] < 200 and 3 < s["e"] < 7          # chiral scale, Skyrme coupling
    assert 700 < s["M_cl"] < 1000                           # classical soliton mass ~ baryon scale
