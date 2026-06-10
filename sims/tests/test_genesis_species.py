"""Species-resolved genesis: batched Skyrme sectors, masses from soliton dynamics, one transition."""
import numpy as np
from cartasis_sims import genesis_species as gs


def test_species_masses_rise_with_coupling():
    # one knot per sector, relaxed: heavier coupling -> heavier knot, ~ sqrt(kappa) (Skyrme scaling)
    species = [("a", 12.0, 0.015), ("b", 24.0, 0.010), ("c", 48.0, 0.006)]
    r = gs.relax_solitons(species, L=28, steps=300, dt=0.004)
    assert list(r["mass_rel"]) == sorted(r["mass_rel"])     # ordered by coupling
    assert r["mass_rel"][-1] > 1.4                          # a real spread, not flat
    assert np.all(np.abs(r["B"]) > 0.7)                     # each knot holds its baryon number


def test_species_condense_together_batched():
    # several sectors cool through ONE transition and each condenses its species (batched in one pass)
    species = [("light", 18.0, 0.012), ("heavy", 40.0, 0.007)]
    res = gs.run_species(species, L=24, steps=500, dt=0.005, seed=1, record=4)
    assert res["energy"].shape[1] == 2                      # both sectors evolved in the batch
    assert np.all(res["energy"][-1] < 0.3 * res["energy"][0])   # both condensed (energy collapses)
    assert np.all(res["content"][-1] > 0.5)                 # both left a relic of knots
