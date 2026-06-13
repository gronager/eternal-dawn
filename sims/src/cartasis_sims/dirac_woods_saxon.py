r"""Dirac in a Woods-Saxon bag: how many radial generations a torsiton bag actually binds.

THE MISSING ANALYTIC LINK. Eq.~6.2 (the nonlinear Dirac equation with the torsion four-fermion
term) in mean field is a Dirac particle of vacuum (constituent) mass m_vac moving in the
self-consistent SCALAR well M(r) -- the chiral-restored bag (M -> 0 in the core where the
condensate melts, M -> m_vac outside). Approximate that well by a Woods-Saxon profile -- the
generic shape of any saturating, surface-thickness-`a` bag of half-density radius R0:

    M(r) = m_vac / (1 + exp(-(r - R0)/a))            # ~0 inside the bag, m_vac outside

The bound radial levels of the lowest channel (kappa=-1) in THIS finite well are the conjectured
generations. The lepton-span fit *assumed* there are three and tuned s_T to their ratios; this
module instead ASKS the well: how many levels does a bag of a given sharpness s_T = R0/r0 actually
bind, over what range of s_T is the answer exactly three (natural vs fine-tuned), and is that s_T
the lattice value (run/09: ~0.36; run/10: ~0.86) -- or nowhere near it?

It does NOT re-postulate the count: n_bound is an OUTPUT of the Dirac solver
(dirac_soliton.dirac_levels), not an input. This is the honest test the hierarchy story owes.

Units: r0 = 1 (lengths in Sommer units), so R0 = s_T directly; m_vac is the constituent mass in
1/r0 (QCD-like ~0.9). m_vac*r0 is the ONE substrate combination -- vary it to see the sensitivity.
"""
from __future__ import annotations

import numpy as np

from . import dirac_soliton as ds

_trapz = getattr(np, "trapezoid", None) or np.trapz

# observed charged-lepton masses (MeV) for the ratio comparison
LEPTONS = np.array([0.510999, 105.658, 1776.86])     # e, mu, tau
LEPTON_SPAN = LEPTONS[-1] / LEPTONS[0]               # ~3477


def ws_mass(r, m_vac, R0, a):
    """The Woods-Saxon scalar bag M(r): ~0 in the core (r<R0), m_vac outside, surface width a."""
    return m_vac / (1.0 + np.exp(-(r - R0) / a))


def bag_levels(s_T, m_vac=0.9, a_frac=0.15, kappa=-1, R_mult=12.0, N=900, nmax=12):
    """Bound radial Dirac levels of a Woods-Saxon bag of sharpness s_T = R0/r0 (r0=1).

    Returns dict: R0, a, energies E_n (ascending, in units of m_vac... actually 1/r0), n_bound,
    the (G,F) wavefunctions, the well M(r) and grid r. a = a_frac * R0 (surface as a fraction of
    the radius -- the bag's relative sharpness)."""
    R0 = float(s_T)                                   # r0 = 1
    a = max(a_frac * R0, 1e-3)
    R = max(10.0, R_mult * R0)
    r = np.linspace(R / N, R, N)
    M = ws_mass(r, m_vac, R0, a)
    V0 = np.zeros_like(r)
    lv, Mvac = ds.dirac_levels(M, V0, r, kappa=kappa, n_levels=nmax)
    E = np.array([e for e, _, _ in lv])
    return {"R0": R0, "a": a, "m_vac": m_vac, "E": E, "n_bound": len(lv),
            "levels": lv, "M": M, "r": r}


def levels_vs_sT(sT_grid, m_vac=0.9, a_frac=0.15, **kw):
    """n_bound as a function of bag sharpness s_T -- the level-counting 'phase diagram'."""
    return np.array([bag_levels(s, m_vac=m_vac, a_frac=a_frac, **kw)["n_bound"] for s in sT_grid])


def sT_window_for(n_target=3, m_vac=0.9, a_frac=0.15, sT_lo=0.1, sT_hi=6.0, ns=120, **kw):
    """The s_T interval over which the bag binds EXACTLY n_target levels (its width = how
    fine-tuned 'n_target generations' is). Returns (s_T_low, s_T_high, width) or None."""
    grid = np.linspace(sT_lo, sT_hi, ns)
    nb = levels_vs_sT(grid, m_vac=m_vac, a_frac=a_frac, **kw)
    hit = grid[nb == n_target]
    if hit.size == 0:
        return None
    return float(hit.min()), float(hit.max()), float(hit.max() - hit.min())


def configmass_ratios(res):
    """Configurational masses m_n = integral (G_n^2 + F_n^2) M(r) dr (the field energy each rung
    carries in the mass-giving region, the structure of Eq. configmass), normalised to the lightest.
    Returns the ascending ratios m_n/m_0 -- the predicted generation hierarchy from this bag."""
    M, r = res["M"], res["r"]
    m = []
    for e, G, F in res["levels"]:
        m.append(_trapz((G**2 + F**2) * M, r))
    m = np.array(m)
    if m.size == 0:
        return m
    order = np.argsort([e for e, _, _ in res["levels"]])
    m = m[order]
    return m / m[0]


def overlap_ratios(res, source_size=0.25):
    """The 'Yukawa-as-overlap' generation masses against a SHARP central source
    chi(r)=exp(-(r/source_size)^2): m_n ~ |<chi|G_n>|^2. A sharp source makes higher rungs (which
    extend outward) overlap exponentially less -> a steep hierarchy. BUT the sharpness is an EXTRA
    assumption (the would-be-Higgs core), not the physical configmass, which uses the broad M(r)."""
    r = res["r"]
    chi = np.exp(-(r / source_size) ** 2)
    o, Es = [], []
    for e, G, F in res["levels"]:
        norm = _trapz(G**2 + F**2, r)
        o.append(abs(_trapz(G * chi, r)) / np.sqrt(max(norm, 1e-300)))
        Es.append(e)
    o = np.array(o)
    if o.size == 0:
        return o
    m = (o[np.argsort(Es)]) ** 2
    return m / np.where(m[0] > 0, m[0], 1.0)


def report(m_vac=0.9, a_frac=0.15, lattice_sT=(0.36, 0.86)):
    """Print the honest verdict: the level-count phase diagram, the s_T that yields 3 generations,
    its width (natural vs fine-tuned), the predicted mass ratios there, and where the lattice sits."""
    print(f"Dirac in a Woods-Saxon bag  (m_vac*r0 = {m_vac:.2f}, surface a = {a_frac:.2f} R0)\n")
    print("  s_T = R0/r0   n_bound   lowest E_n / m_vac")
    for s in (0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0):
        res = bag_levels(s, m_vac=m_vac, a_frac=a_frac)
        e = res["E"] / m_vac
        estr = " ".join(f"{x:.2f}" for x in e[:4])
        print(f"    {s:>6.2f}      {res['n_bound']:>3d}      {estr}")

    win = sT_window_for(3, m_vac=m_vac, a_frac=a_frac)
    print()
    if win:
        lo, hi, w = win
        nat = "a BROAD, natural band" if w > 0.5 else "a NARROW, fine-tuned sliver"
        print(f"  exactly 3 bound levels for s_T in [{lo:.2f}, {hi:.2f}]  (width {w:.2f}) -- {nat}")
        mid = 0.5 * (lo + hi)
        res3 = bag_levels(mid, m_vac=m_vac, a_frac=a_frac)
        ratios = configmass_ratios(res3)
        sharp = overlap_ratios(res3, source_size=0.25)
        if ratios.size >= 3:
            print(f"  at s_T={mid:.2f}: configmass (broad M) span {ratios[2]/ratios[0]:.1f}; "
                  f"sharp-source overlap span {sharp[2]/sharp[0]:.1f}  "
                  f"-- both vs observed lepton span {LEPTON_SPAN:.0f}")
            print(f"    => the bag-level mechanism gives the right SHAPE but ~O(1-10) magnitude, "
                  f"not ~{LEPTON_SPAN:.0f}: the hierarchy is not a derived output here.")
    else:
        print("  NO s_T in [0.1, 6] binds exactly 3 levels at this (m_vac, a).")
    lo_l, hi_l = lattice_sT
    nb_l = [bag_levels(s, m_vac=m_vac, a_frac=a_frac)["n_bound"] for s in (lo_l, hi_l)]
    print(f"\n  the lattice bag: s_T(09)={lo_l} -> {nb_l[0]} level(s);  "
          f"s_T(10)={hi_l} -> {nb_l[1]} level(s)")
    verdict = ("the lattice bag is large enough for 3" if min(nb_l) >= 3 else
               "the lattice bag binds FEWER than 3 -- the generations need a deeper well "
               "(stronger torsion coupling) or a bigger bag than measured")
    print(f"  VERDICT: {verdict}.")


def main():
    report()


if __name__ == "__main__":
    main()
