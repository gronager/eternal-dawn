r"""Horizon, step 3d: the unified genesis -- ALL observable matter from ONE cooling transition.

The capstone synthesis. 3a (confinement), 3b (baryons), 3c (leptons) become one timeline: as the
chiral condensate v(T) switches on through a single transition, every species acquires mass AT ONCE
and the matter content of the universe condenses out --
  * the COLOURED sector -> quarks, confined (3a) into the baryon octet/decuplet (3b);
  * the COLOURLESS sector -> the lepton generations e,mu,tau and the sub-eV neutrinos (3c).
One condensate, two organising principles, the whole Standard-Model matter spectrum in one quench.

PARAMETRISED FOR THE LATTICE. The single owed input is the scale and the bag sharpness. `genesis()`
takes `scale_MeV` (= Lambda, the condensation scale) and `s_T` (the generation rung) -- defaults are
the calibrated/illustrative values so it runs NOW; pass the lattice's Lambda and s_T when the L4
ensemble lands and the *same* call returns the absolute-MeV prediction. This is the re-run hook.

HONEST: baryons ~30% (Skyrme + rigid rotor); leptons from the soliton ladder + one s_T; neutrino
masses illustrative; the abundance is the Kibble--Zurek FORMATION density (how much matter a cooling
bounce makes -- the full baryon asymmetry is the genesis-chapter physics, not this). The *structure*
-- one transition, all species, the spectrum -- is the content; the *numbers* are the lattice's word."""
from __future__ import annotations

import numpy as np

from . import genesis_su3 as su3
from . import genesis_leptons as lep


def full_matter(scale_MeV=None, s_T=None):
    """The complete observable matter spectrum from the soliton: the baryon octet+decuplet (3b) and the
    lepton generations + neutrinos (3c). `scale_MeV` rescales the whole spectrum to the lattice Lambda
    (default: the N/Delta/Omega-calibrated MeV); `s_T` re-rungs the charged leptons (default: fitted).
    Returns {species: (mass_MeV, sector)} with sector in {octet, decuplet, lepton, neutrino}."""
    bar = su3.baryon_spectrum()["mass"]                    # calibrated MeV
    charged_names = ["e", "mu", "tau"]
    if s_T is not None:                                    # re-rung the leptons at the lattice s_T
        from . import fermion_masses as fm
        lad = fm._ascending_ladder(s_T, well=fm.WELL)      # ascending ladder, tau = the largest rung
        charged = lad / lad[-1] * 1776.86                  # anchor the tau at its observed mass (MeV)
    else:
        charged = lep.charged_lepton_tower()["mass"]
    nu = lep.neutrino_tower()

    out = {}
    for b in ("N", "Lambda", "Sigma", "Xi"):
        out[b] = (float(bar[b]), "octet")
    for b in ("Delta", "Sigma*", "Xi*", "Omega"):
        out[b] = (float(bar[b]), "decuplet")
    for nm, m in zip(charged_names, np.asarray(charged, float)):
        out[nm] = (float(m), "lepton")
    for nm, m in zip(nu["names"], nu["mass"]):
        out[nm] = (float(m), "neutrino")

    if scale_MeV is not None:                              # rescale to the lattice condensation scale
        ref = out["N"][0]                                  # anchor on the nucleon (= the soliton scale)
        f = scale_MeV / ref
        out = {k: (m * f, s) for k, (m, s) in out.items()}
    return out


def condensate_curve(temps, Tc=0.5, width=0.05):
    """The chiral order parameter v(T): 0 hot (symmetric) -> 1 cold (condensed vacuum), a tanh ramp."""
    return 0.5 * (1.0 + np.tanh((Tc - np.asarray(temps, float)) / width))


def kibble_zurek_abundance(cooling_rate, sector):
    """Schematic Kibble--Zurek FORMATION density per sector: faster quench -> more knots (defect
    density ~ rate^nu). Coloured (baryon) and colourless (lepton) knots both scale; the relative number
    is the genesis matter content. Returns a relative abundance (illustrative -- the asymmetry is owed)."""
    nu = 0.5                                               # mean-field KZ exponent (illustrative)
    base = cooling_rate ** nu
    weight = {"octet": 1.0, "decuplet": 0.5, "lepton": 1.0, "neutrino": 1.0}.get(sector, 1.0)
    return base * weight


def genesis(scale_MeV=None, s_T=None, Tc=0.5, width=0.05, cooling_rate=1.0, n_T=240):
    """The unified genesis: all observable matter condensing from one transition. Returns the cooling
    timeline (temps, v), every species' mass switching on as v(T) grows, the final spectrum, and the
    Kibble--Zurek formation abundances. PASS scale_MeV (=Lambda) and s_T from the lattice to re-run for
    the absolute-MeV prediction; defaults are the calibrated/illustrative values (runs now)."""
    temps = np.linspace(1.0, 0.0, n_T)
    v = condensate_curve(temps, Tc=Tc, width=width)
    matter = full_matter(scale_MeV=scale_MeV, s_T=s_T)
    masses = {sp: v * m for sp, (m, _) in matter.items()}  # m(T) = v(T) * m_final, switching on together
    abundance = {sp: kibble_zurek_abundance(cooling_rate, sec) for sp, (_, sec) in matter.items()}
    return {"temps": temps, "v": v, "masses": masses, "matter": matter,
            "abundance": abundance, "Tc": Tc,
            "scale_MeV": scale_MeV, "s_T": s_T, "cooling_rate": cooling_rate}
