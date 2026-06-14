r"""Hopfions: S^2-valued textures with nonzero Hopf charge (self-linking) = the Route-1 generation
candidate. Stage 1 here: CONSTRUCT charge-Q Hopfion fields, VERIFY the Hopf invariant numerically,
and measure the Faddeev-Skyrme energy ladder E(Q) -- the candidate generation mass ladder.

A generation is conjectured to be a distinct linking number Q = 0,1,2 of the neutral order-parameter
field n_hat(x) in S^2 (pi_3(S^2)=Z). Distinct Q are distinct topological sectors -> separately
protected (no inter-generation transition, mu->e gamma forbidden) AND an ordered ladder (this module).

Construction: inverse stereographic R^3 -> S^3 in C^2,
    Z1 = 2(x+iy)/(1+r^2),   Z2 = (-2z + i(1-r^2))/(1+r^2),   |Z1|^2+|Z2|^2 = 1,
then the Hopf base coordinate w = Z1^a / Z2^b (Hopf charge Q = a*b), and n_hat from w by stereographic
    n_hat = (2 Re w, 2 Im w, |w|^2 - 1)/(|w|^2 + 1)   -> (0,0,1) at infinity (the vacuum).

Hopf charge (numeric): f_jk = n_hat . (d_j n_hat x d_k n_hat) is the pullback area 2-form; B_i =
1/2 eps_ijk f_jk satisfies div B = 0, so B = curl A (solve A in Fourier, Coulomb gauge); then
    Q = (1/(16 pi^2)) integral A . B d^3x      (calibrated so the elementary map gives Q=1).
Energy (Faddeev-Skyrme):  E = integral [ (d_i n_hat).(d_i n_hat) + 1/2 f_jk f_jk ] d^3x.
The Vakulenko-Kapitanskii bound E >= c |Q|^{3/4} sets the expected ladder.
"""
from __future__ import annotations

import numpy as np

_pi = np.pi


def grid(L=4.0, N=64):
    """Periodic box [-L,L)^3, N^3 points. Returns (x,y,z) meshes, spacing h, and k-meshes for FFT."""
    ax = np.linspace(-L, L, N, endpoint=False)
    h = ax[1] - ax[0]
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    k1 = 2 * _pi * np.fft.fftfreq(N, d=h)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    return X, Y, Z, h, (KX, KY, KZ)


def hopfion_field(X, Y, Z, a=1, b=1):
    """n_hat(x) for a Hopfion of charge Q=a*b. Inverse-stereographic R^3->S^3 (Z1,Z2), then the Hopf
    map of (Z1^a, Z2^b) in the DIVISION-FREE bilinear form (regular everywhere, no w=Z1/Z2 blow-up):
        n_hat = (2 Re(conj(z1) z2), 2 Im(conj(z1) z2), |z1|^2-|z2|^2)/(|z1|^2+|z2|^2),  z1=Z1^a, z2=Z2^b."""
    r2 = X**2 + Y**2 + Z**2
    den = 1.0 + r2
    Z1 = (2 * (X + 1j * Y)) / den
    Z2 = (-2 * Z + 1j * (1.0 - r2)) / den
    z1, z2 = Z1**a, Z2**b
    Nsq = np.abs(z1) ** 2 + np.abs(z2) ** 2
    zz = np.conj(z1) * z2
    n1 = 2 * np.real(zz) / Nsq
    n2 = 2 * np.imag(zz) / Nsq
    n3 = (np.abs(z1) ** 2 - np.abs(z2) ** 2) / Nsq
    n = np.stack([n1, n2, n3], axis=0)
    return n / np.linalg.norm(n, axis=0)


def _d(field, K, axis):
    """Spectral derivative d/dx_axis of a real scalar field on the periodic box."""
    return np.real(np.fft.ifftn(1j * K[axis] * np.fft.fftn(field)))


def _fjk(n, K):
    """Antisymmetric f_jk = n.(d_j n x d_k n); returns the three independent comps f_01,f_02,f_12."""
    dn = [[_d(n[a], K, j) for j in range(3)] for a in range(3)]   # dn[a][j] = d_j n_a

    def f(j, k):
        # n . (d_j n x d_k n)
        djn = np.stack([dn[a][j] for a in range(3)], axis=0)
        dkn = np.stack([dn[a][k] for a in range(3)], axis=0)
        cross = np.cross(djn, dkn, axis=0)
        return np.sum(n * cross, axis=0)

    return f(0, 1), f(0, 2), f(1, 2), dn


def hopf_charge(n, K, h):
    """Numeric Hopf invariant. B = (f_12, -f_02, f_01) = curl A; solve A in Coulomb gauge by FFT;
    Q = (1/16 pi^2) integral A.B. (Calibrated to Q=1 for the elementary a=b=1 map.)"""
    f01, f02, f12, _ = _fjk(n, K)
    B = [f12, -f02, f01]                                  # B_i = 1/2 eps_ijk f_jk
    Bk = [np.fft.fftn(Bi) for Bi in B]
    KX, KY, KZ = K
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0
    # A = curl^{-1} B:  in Fourier, B = i k x A (Coulomb gauge k.A=0) -> A = -i (k x B)/k^2
    kxB = [KY * Bk[2] - KZ * Bk[1], KZ * Bk[0] - KX * Bk[2], KX * Bk[1] - KY * Bk[0]]
    Ak = [-1j * c / k2 for c in kxB]
    for c in Ak:
        c[0, 0, 0] = 0.0
    A = [np.real(np.fft.ifftn(c)) for c in Ak]
    AB = sum(A[i] * B[i] for i in range(3))
    Q = (1.0 / (16 * _pi**2)) * np.sum(AB) * h**3
    return float(Q)


def faddeev_energy(n, K, h):
    """Faddeev-Skyrme energy E = integral [ (d n)^2 + 1/2 f_jk f_jk ] d^3x. (E2 = gradient, E4 = Skyrme.)"""
    f01, f02, f12, dn = _fjk(n, K)
    E2 = sum(np.sum(dn[a][j] ** 2) for a in range(3) for j in range(3)) * h**3
    E4 = np.sum(f01**2 + f02**2 + f12**2) * h**3          # 1/2 * f_jk f_jk over j<k pairs *2 = sum sq
    return float(E2), float(E4), float(E2 + E4)


def report(L=4.0, N=64):
    print(f"Hopfion ladder (box L={L}, N={N}^3):  Q = linking number = generation candidate\n")
    print(f"  {'(a,b)':>7} {'Q_numeric':>10} {'E2':>9} {'E4':>9} {'E_total':>9} {'E/Q^{3/4}':>10}")
    X, Y, Z, h, K = grid(L, N)
    base = None
    for (a, b) in [(1, 1), (1, 2), (1, 3), (2, 1)]:
        n = hopfion_field(X, Y, Z, a, b)
        Q = hopf_charge(n, K, h)
        E2, E4, E = faddeev_energy(n, K, h)
        qabs = max(abs(Q), 1e-9)
        print(f"  {f'({a},{b})':>7} {Q:>10.3f} {E2:>9.2f} {E4:>9.2f} {E:>9.2f} {E/qabs**0.75:>10.2f}")
        if base is None:
            base = E
    print("\n  Q should come out ~ a*b; E should ladder ~ Q^{3/4} (Vakulenko-Kapitanskii).")
    print("  The mass ladder E(Q) is the candidate generation hierarchy (Stage 1: field energies only;")
    print("  the Dirac zero modes in each background are Stage 2).")


if __name__ == "__main__":
    report()
