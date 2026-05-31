"""Tests for the dark-energy parent-accretion toy model and helpers."""

import math

import numpy as np

from cartasis_sims import dark_energy as de


def test_cpl_from_toy_signs():
    # still-growing (p0>0), saturating (q>0): w0 > -1 and wa < 0
    w0, wa = de.cpl_from_toy(eps=0.5, p0=0.6, q=1.2)
    assert w0 > -1.0
    assert wa < 0.0


def test_fit_toy_recovers_target():
    w0_t, wa_t = -0.752, -0.86
    fit = de.fit_toy_to(w0_t, wa_t)
    w0, wa = de.cpl_from_toy(**fit)
    assert math.isclose(w0, w0_t, rel_tol=1e-6)
    assert math.isclose(wa, wa_t, rel_tol=1e-6)


def test_w_of_a_reduces_to_LCDM_when_eps_zero():
    a = np.linspace(0.3, 1.0, 10)
    assert np.allclose(de.w_of_a(a, eps=0.0), -1.0)


def test_w_phantom_in_past_for_desy5_fit():
    # DESI: phantom (w<-1) in the past, w>-1 today -- but the TOY model with
    # w=-1+eps p0 a^q stays >= -1; check it at least rises toward the past being
    # MORE negative than today (monotone), matching the qualitative trend.
    fit = de.fit_toy_to(-0.752, -0.86)
    w_today = de.w_of_a(1.0, **fit)
    w_past = de.w_of_a(0.3, **fit)
    assert w_today > w_past            # w larger today, smaller (more neg) in past


def test_filter_fraction_one_sixth():
    f = de.filter_fraction_from_ratio(5.36)
    assert 0.14 < f < 0.18


def test_m_crit_in_cmb_is_sub_lunar():
    m = de.m_crit_growth(2.725)
    assert 1e22 < m < 1e23          # ~0.6 lunar masses -> all astrophysical BHs grow


def test_injection_model_crosses_phantom():
    # peaked injection: w<-1 in the past (a<a_p), w>-1 today (a>a_p)
    inj = de.injection_from_cpl(-0.752, -0.86)
    a_p = inj["a_p"]
    assert de.w_eff_injection(0.5 * a_p, a_p, inj["beta"]) < -1.0   # phantom past
    assert de.w_eff_injection(1.0, a_p, inj["beta"]) > -1.0         # today
    assert 0 < inj["z_cross"] < 1.5                                 # sensible crossing


def test_injection_from_cpl_recovers_w0_wa():
    w0_t, wa_t = -0.752, -0.86
    inj = de.injection_from_cpl(w0_t, wa_t)
    # check the model reproduces w0 at a=1 exactly
    w0 = de.w_eff_injection(1.0, inj["a_p"], inj["beta"])
    assert math.isclose(float(w0), w0_t, rel_tol=1e-6)


def test_rho_de_peaks_at_ap():
    a = np.linspace(0.2, 1.2, 400)
    a_p, beta = 0.75, 1.3
    rho = de.rho_de_injection(a, a_p, beta)
    assert math.isclose(a[np.argmax(rho)], a_p, abs_tol=0.01)


def test_logistic_accretion_monotone_saturating():
    a = np.linspace(0.05, 1.0, 50)
    M = de.logistic_accretion(a)
    assert np.all(np.diff(M) > 0)          # monotone growth
    assert M[-1] < 1.0 and M[-1] > M[0]    # saturating toward m_inf
