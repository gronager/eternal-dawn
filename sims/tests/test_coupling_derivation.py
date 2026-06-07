"""SymPy derivation of the torsiton couplings -- the relations, verified symbolically."""
import sympy as sp

from cartasis_sims import coupling_derivation as cd


def test_coupling_ratios_symbolic():
    g = sp.symbols("g", positive=True)
    rel = cd.coupling_relations()
    assert sp.simplify(rel["lambda"] - 2 * g**2) == 0          # lambda = 2 g^2
    assert sp.simplify(rel["g_v"] - sp.sqrt(3) * g) == 0       # g_v = sqrt(3) g
    assert rel["GV_over_GS"] == 2                               # the Fierz ratio, recovered


def test_loop_integrals_have_the_right_uv_behaviour():
    # I0 ~ Lambda^2 (quadratic), I1 ~ log(Lambda) -- the standard NJL divergences
    M, Lam2, I0, I1 = cd.loop_integrals()
    big = {M: 1}
    # leading large-Lambda^2: I0 grows like Lambda2, I1 like log(Lambda2)
    I0_lead = sp.limit(I0.subs(big) / Lam2, Lam2, sp.oo)
    assert I0_lead == 1 / (16 * sp.pi**2)                       # I0 ~ Lambda^2/16pi^2
    assert sp.limit(I1.subs(big), Lam2, sp.oo) == sp.oo         # I1 diverges (logarithmically)


def test_sea_condensate_is_negative_and_grows_with_mass():
    # <qbar q> = -4 N_c M I0 < 0, and |<qbar q>| increases with M at fixed cutoff
    M, Lam2, I0, I1 = cd.loop_integrals()
    qq = cd.sea_condensate()
    num = qq.subs({sp.symbols("N_c", positive=True): 3, Lam2: 9})
    v1 = float(num.subs(M, 0.5)); v2 = float(num.subs(M, 1.0))
    assert v1 < 0 and v2 < 0 and abs(v2) > abs(v1)
