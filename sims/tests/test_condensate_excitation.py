"""Generations as condensate excitations: ceiling-free span, near-orthogonality node, Koide-compatible."""
import numpy as np

from cartasis_sims import condensate_excitation as ce


def test_span_is_ceiling_free_near_orthogonality():
    # the span diverges as the fermion approaches orthogonality with the excited modes (a/b -> 1):
    # the magnitude mechanism the bound-state towers (capped ~2-15) lacked
    assert ce.span(0.99) > 1000                  # easily reaches/exceeds the lepton span 3477
    assert ce.span(0.99) > 100 * ce.span(0.85)   # blows up toward orthogonality


def test_lightest_is_most_excited_mode():
    # the smallest |overlap| is the highest condensate excitation (most cancellation) = light gen
    m = np.abs(ce.overlap_masses(0.85))
    assert np.argmin(m) == len(m) - 1            # n=2 mode is the lightest


def test_koide_crossed_by_overlap_masses():
    # Q(a/b) is continuous and brackets 2/3 -> the overlap masses are Koide-compatible at some a/b
    assert ce.koide_Q(0.95) < 2 / 3 < ce.koide_Q(0.97)


def test_criticality_drives_span_up():
    # approaching marginal binding (the chiral boundary) grows the span and the binding margin:
    # the criticality driver is real (span larger at the more-marginal coupling)
    sc = ce.criticality_scan(g_values=[4.0, 2.7])
    deep, marginal = sc[0], sc[1]
    assert marginal["margin"] > deep["margin"]      # E/Mvac -> 1 (more marginal)
    assert marginal["span"] > 2 * deep["span"]      # span grows sharply toward criticality
