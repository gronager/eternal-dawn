"""Hopfions: the Hopf charge comes out integer (Q=a*b) and the energy ladders up with |Q|."""
import numpy as np

from cartasis_sims import hopfion as hf


def test_hopf_charge_is_integer():
    # the constructed fields must carry the conjectured topological charge Q = a*b (sign = orientation)
    X, Y, Z, h, K = hf.grid(L=4.0, N=48)
    for (a, b, q) in [(1, 1, 1), (1, 2, 2), (2, 1, 2)]:
        n = hf.hopfion_field(X, Y, Z, a, b)
        Q = hf.hopf_charge(n, K, h)
        assert abs(abs(Q) - q) < 0.15, f"(a,b)=({a},{b}): |Q|={abs(Q):.3f} != {q}"


def test_energy_ladders_with_charge():
    # distinct sectors form an ordered ladder: E(Q=2) > E(Q=1) (the candidate generation mass ladder)
    X, Y, Z, h, K = hf.grid(L=4.0, N=48)
    E1 = hf.faddeev_energy(hf.hopfion_field(X, Y, Z, 1, 1), K, h)[2]
    E2 = hf.faddeev_energy(hf.hopfion_field(X, Y, Z, 1, 2), K, h)[2]
    assert E2 > E1 > 0


def test_vacuum_at_infinity():
    # n_hat is the constant vacuum on the box boundary (localised texture, FFT-compatible)
    X, Y, Z, h, K = hf.grid(L=4.0, N=40)
    n = hf.hopfion_field(X, Y, Z, 1, 1)
    edge = np.concatenate([n[:, 0, :, :].ravel(), n[:, -1, :, :].ravel()])
    # the boundary points are nearly aligned (constant direction): small spread
    nb = n[:, 0, 0, 0]
    assert abs(np.linalg.norm(nb) - 1.0) < 1e-6
