r"""The OGU as a finite-lived recursive timeline: saturation, evaporation, depth.

This resolves the apparent tension between "OGUs grow toward Nariai" and "OGUs are
immortal (~1e142 s)". Both hold, in sequence, and the timescale hierarchy makes the
recursion well-defined.

1. Growth and SATURATION. An OGU grows fast-then-slow (settler regime, void_scale.py)
   toward the Nariai mass in ~1 Hubble time. But growth does not continue forever: a
   black hole's Hawking temperature T_BH = hbar c^3/(8 pi G M k_B) FALLS as it grows,
   and the de Sitter void it sits in has a floor temperature T_dS = hbar H/(2 pi k_B).
   They equalise at

       M_eq = c^3/(4 G H_Lambda) ~ 5.6e52 kg  ~  M_Nariai ~ a Hubble mass,

   where the OGU is in thermal equilibrium with its own cosmological horizon: no net
   flow, growth stops. So an OGU SATURATES at ~a Hubble mass (our mass) -- it does not
   keep growing without bound; dark energy caps it (cf. Nariai, void_scale.py).

2. EVAPORATION. At saturation T_BH ~ T_dS, so the OGU is in (unstable) equilibrium and
   then very slowly evaporates on the Hawking time t_evap ~ 5120 pi G^2 M^3/(hbar c^4)
   ~ 1e142 s. It evaporates LONG before it could ever touch a neighbour (separation
   ~ e^{I/4} horizons, utterly causally disconnected), so the membrane is finite-lived
   but astronomically long-lived: it persists ~1e142 s, then returns to the void.

3. RECURSION DEPTH per lifetime. Meanwhile the INTERIOR runs its own timeline: it
   reaches heat death and seeds its own sub-OGUs every recursion "tick" ~ the interior
   black-hole-evaporation era ~1e108 s (1e100 yr). The host OGU lives ~1e142 s, so it
   hosts

       N_ticks ~ t_evap / t_tick ~ 1e34

   full recursive interior cycles before its own membrane evaporates. So the answer to
   "recursive timelines forever, or pulsating in the flat void?" is BOTH, cleanly: each
   membrane is FINITE-lived (it evaporates), so no single branch nests forever; but each
   hosts ~1e34 interior cycles, and new OGUs nucleate without end both in the flat void
   AND inside every heat-dead interior. The tree is eternal and ever-expanding into the
   void even though every node is mortal -- a fractal of finite, pulsating timelines,
   not an actually-infinite nesting. The void is the one truly eternal, flat substrate;
   everything in it is a mortal, recursive bubble.
"""

from __future__ import annotations

import numpy as np

from . import constants as k

T_INTERIOR_TICK = 1.0e108     # s; interior heat-death / BH-evaporation era (~1e100 yr)


def hubble_lambda() -> float:
    return k.H0 * np.sqrt(k.Omega_Lambda)


def hawking_temperature(M: float) -> float:
    """T_BH = hbar c^3/(8 pi G M k_B) [K]; falls as the hole grows."""
    return k.hbar * k.c**3 / (8.0 * np.pi * k.G * M * k.kB)


def de_sitter_temperature() -> float:
    """T_dS = hbar H_Lambda/(2 pi k_B) [K]: the void's floor temperature."""
    return k.hbar * hubble_lambda() / (2.0 * np.pi * k.kB)


def equilibrium_mass() -> float:
    """M_eq where T_BH = T_dS: c^3/(4 G H_Lambda) [kg]. The OGU saturates here,
    in thermal equilibrium with its own cosmological horizon (~a Hubble mass)."""
    return k.c**3 / (4.0 * k.G * hubble_lambda())


def saturates(M: float) -> bool:
    """True once the OGU has reached thermal equilibrium with the void (stops growing)."""
    return M >= equilibrium_mass()


def evaporation_time(M: float) -> float:
    """Hawking evaporation time 5120 pi G^2 M^3/(hbar c^4) [s]."""
    return 5120.0 * np.pi * k.G**2 * M**3 / (k.hbar * k.c**4)


def membrane_lifetime() -> float:
    """How long an OGU's Cartasis membrane persists: t_evap at the saturation mass [s]."""
    return evaporation_time(equilibrium_mass())


def recursion_cycles_per_lifetime(t_tick: float = T_INTERIOR_TICK) -> float:
    """Number of full interior recursive cycles an OGU hosts before its membrane
    evaporates: t_evap(M_eq)/t_tick ~ 1e34."""
    return membrane_lifetime() / t_tick


def touches_neighbour_before_evaporating(separation_horizons: float) -> bool:
    """An OGU would have to grow across `separation_horizons` cosmological horizons to
    reach a neighbour. In the dilute regime separation ~ e^{I/4} >> 1, and growth caps
    at ~1 horizon (Nariai), so it NEVER touches a neighbour -- it evaporates first."""
    return separation_horizons <= 1.0
