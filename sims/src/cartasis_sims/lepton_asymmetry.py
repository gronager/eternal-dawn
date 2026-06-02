r"""The neutrino sector's "unmatched number": the relic lepton asymmetry.

The baryon-to-photon ratio is one of the best-measured numbers in cosmology,
eta = n_B/n_gamma ~ 6.1e-10 (Ch8: doubly anchored by the CMB peaks and BBN
deuterium). The *lepton* asymmetry of the universe -- carried mostly by the
cosmic neutrino background as a neutrino/antineutrino excess -- is, by contrast,
almost completely unmeasured. That gap is the test.

Two facts collide:

1. SCT's prediction. The seed asymmetry is inherited, then processed at the
   bounce by electroweak sphalerons, which Ch4 shows are *active* there
   (T ~ 1e7 GeV >> the ~130 GeV sphaleron freeze-out). Sphalerons in equilibrium
   conserve B-L but wash B and L into a fixed ratio. For the Standard Model with
   n_gen=3 generations and n_H=1 Higgs doublet,

       B =  (8 n_gen + 4 n_H)/(22 n_gen + 13 n_H) (B-L) = 28/79 (B-L),
       L =  B - (B-L)                              = -51/79 (B-L),

   so |L/B| = 51/28 ~ 1.82: the lepton asymmetry is the *same order* as the
   baryon asymmetry, not enhanced. The relic neutrino degeneracy parameter is
   therefore xi_nu ~ eta ~ 1e-9 -- utterly indistinguishable from zero. SCT does
   NOT have a large primordial lepton asymmetry to spend.

2. The observational ceiling. BBN + CMB bound the common neutrino degeneracy
   only to |xi_nu| <~ 0.04 (e.g. Oldengott & Schwarz 2017; Planck 2018). That
   ceiling permits a neutrino/antineutrino asymmetry up to ~1e7 times the baryon
   asymmetry -- a vast unmatched headroom in which a lepton-only asymmetry could
   hide.

The discriminating prediction: SCT lives near the floor (xi_nu ~ eta), not the
ceiling. A *measured* large primordial lepton asymmetry (xi_nu >~ 1e-2) would
require a lepton-number source acting AFTER sphaleron freeze-out -- which the
inherit-then-sphaleron-process picture does not supply -- and would disfavour
SCT's account of the asymmetry. A null result (xi_nu consistent with ~eta) is a
weak confirmation. PTOLEMY-class C-nu-B experiments and improved BBN+CMB
analyses are the route.

Reproducible companion to chapters/08-observational-tests.tex,
"The neutrino sector's unmatched number".
"""

from __future__ import annotations

from scipy.special import zeta

ETA_OBS = 6.1e-10           # baryon-to-photon ratio n_B/n_gamma (doubly anchored)
XI_CEILING = 0.04           # |xi_nu| ceiling from BBN+CMB (Oldengott-Schwarz 2017)
TNU_TGAMMA_CUBED = 4.0 / 11.0   # (T_nu/T_gamma)^3 after e+e- annihilation


def sphaleron_BL_coefficients(n_gen: int = 3, n_H: int = 1) -> tuple[float, float]:
    """(c_B, c_L) such that B = c_B (B-L), L = c_L (B-L) for sphalerons in
    equilibrium. Standard Model: c_B = 28/79, c_L = -51/79 (n_gen=3, n_H=1)."""
    num_B = 8 * n_gen + 4 * n_H
    den = 22 * n_gen + 13 * n_H
    c_B = num_B / den
    c_L = c_B - 1.0           # L = B - (B-L)
    return c_B, c_L


def lepton_to_baryon_ratio(n_gen: int = 3, n_H: int = 1) -> float:
    """|L/B| left by sphaleron equilibrium. SM: 51/28 ~ 1.82 -- order unity, NOT
    an enhancement. The lepton asymmetry is comparable to the baryon asymmetry."""
    c_B, c_L = sphaleron_BL_coefficients(n_gen, n_H)
    return abs(c_L / c_B)


def degeneracy_from_lepton_asymmetry(eta_L: float) -> float:
    """Invert the small-xi relation between the per-species neutrino degeneracy
    parameter xi = mu_nu/T_nu and the neutrino lepton-to-photon ratio

        eta_L,nu = (n_nu - n_nubar)/n_gamma
                 ~ (pi^2 / (12 zeta(3))) (T_nu/T_gamma)^3 xi          (small xi),

    returning xi for a given eta_L,nu. For eta_L ~ eta_B this gives xi ~ 1e-9."""
    import math
    prefac = (math.pi**2 / (12.0 * zeta(3))) * TNU_TGAMMA_CUBED
    return eta_L / prefac


def predicted_neutrino_degeneracy(eta_B: float = ETA_OBS,
                                  n_gen: int = 3, n_H: int = 1) -> float:
    """SCT's predicted relic neutrino degeneracy: the sphaleron-locked lepton
    asymmetry eta_L = |L/B| eta_B, converted to xi. ~1e-9, i.e. ~ eta_B."""
    eta_L = lepton_to_baryon_ratio(n_gen, n_H) * eta_B
    return degeneracy_from_lepton_asymmetry(eta_L)


def unmatched_headroom(eta_B: float = ETA_OBS) -> float:
    """How many times larger than SCT's prediction the observational ceiling sits:
    xi_ceiling / xi_predicted. ~1e7 -- the size of the unprobed lepton-asymmetry
    window the framework predicts is empty."""
    return XI_CEILING / predicted_neutrino_degeneracy(eta_B)


def is_discriminating(eta_B: float = ETA_OBS) -> bool:
    """The test discriminates iff SCT's prediction sits far below the current
    ceiling -- so a future detection of a large asymmetry would be decisive."""
    return unmatched_headroom(eta_B) > 1e3
