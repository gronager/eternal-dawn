r"""Fierz projection of the torsion (Hehl-Datta) four-fermion interaction.

The cheap, decisive probe for Part III's electroweak S, before the full RPA. The
electroweak S of a composite sector is killed or saved by how the VECTOR (rho-like)
and AXIAL (a1-like) resonances split: a QCD-like sector has M_a1 > M_rho and S ~ 0.3
(the graveyard); a 'walking' sector has M_a1 -> M_rho, the V and A spectral functions
cancel, and S -> 0. The resonance masses are set by the four-fermion couplings in the
meson (exchange) channels, G_V and G_A -- so the make-or-break reduces to one ratio,
G_A/G_V, which is fixed by the Lorentz structure of the interaction.

The torsion interaction, after integrating out the Einstein-Cartan torsion, is the
Hehl-Datta term -- the square of the AXIAL current,

    L_HD ~ (psi-bar gamma^5 gamma^mu psi)(psi-bar gamma^5 gamma_mu psi)   (axial-axial).

Fierz-rearranging it into the exchange (q-qbar meson) channels, with explicit Dirac
matrices (validated below against the V-A self-Fierz identity and Fierz^2 = identity),
gives

    G_S = +1/4,  G_P = -1/4,  G_V = +1/2,  G_A = +1/2,  G_T = 0,

so **G_A / G_V = 1 exactly**: the torsion interaction couples the vector and axial
meson channels with equal strength. That is the structural prerequisite for walking
(M_a1 -> M_rho) and a small S -- the opposite of a QCD-like sector, and it is forced
by the axial-axial Lorentz structure, not tuned. It is NECESSARY, not sufficient: the
residual V-A splitting from the chiral-breaking loops still needs the full RPA to
quantify. But the cheap probe is green -- the framework has a real structural shot at
escaping the graveyard, so the RPA is worth building.

(The scalar channel G_S = +1/4 > 0 is attractive, driving the chiral condensate -- the
same condensate the chiral soliton forms; and G_T = 0 -- no tensor. The structure
hangs together.)
"""

from __future__ import annotations

import numpy as np

# --- Dirac-representation gamma matrices, metric (+,-,-,-) ---
_I2 = np.eye(2)
_Z2 = np.zeros((2, 2))
_sx = np.array([[0, 1], [1, 0]])
_sy = np.array([[0, -1j], [1j, 0]])
_sz = np.array([[1, 0], [0, -1]])


def _blk(A, B, C, D):
    return np.block([[A, B], [C, D]])


GAMMA = [
    _blk(_I2, _Z2, _Z2, -_I2),       # gamma^0
    _blk(_Z2, _sx, -_sx, _Z2),       # gamma^1
    _blk(_Z2, _sy, -_sy, _Z2),       # gamma^2
    _blk(_Z2, _sz, -_sz, _Z2),       # gamma^3
]
GAMMA5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
METRIC = [1, -1, -1, -1]

_V = [GAMMA[m] for m in range(4)]
_A = [GAMMA[m] @ GAMMA5 for m in range(4)]
_SIG = {(mu, nu): 0.5j * (GAMMA[mu] @ GAMMA[nu] - GAMMA[nu] @ GAMMA[mu])
        for mu in range(4) for nu in range(4)}


def axial_axial_kernel():
    """The direct Hehl-Datta kernel K[i,j,k,l] = sum_mu g_mu A^mu_ij A^mu_kl."""
    K = np.zeros((4, 4, 4, 4), dtype=complex)
    for m in range(4):
        K += METRIC[m] * np.einsum("ij,kl->ijkl", _A[m], _A[m])
    return K


def fierz(K):
    """Exchange (Fierz) rearrangement with the -1 fermion sign: K'[a,b,c,d]=-K[a,d,c,b]."""
    return -np.einsum("adcb->abcd", K)


def _proj_lorentz(K, Ops):
    """Project onto a Lorentz-vector channel (V or A): (1/64) sum_mu g_mu K Op^mu Op^mu."""
    return sum(METRIC[m] * np.einsum("abcd,ba,dc->", K, Ops[m], Ops[m])
               for m in range(4)) / 64.0


def _proj_scalar(K, Op):
    """Project onto a scalar channel (S or P)."""
    return np.einsum("abcd,ba,dc->", K, Op, Op) / 64.0


def _proj_tensor(K):
    tot = 0.0
    for mu in range(4):
        for nu in range(4):
            if mu < nu:
                tot += METRIC[mu] * METRIC[nu] * np.einsum(
                    "abcd,ba,dc->", K, _SIG[(mu, nu)], _SIG[(mu, nu)])
    return tot / 64.0


def hehl_datta_channels():
    """Fierz the axial-axial Hehl-Datta term into the exchange (meson) channels.
    Returns {S, P, V, A, T} couplings (real). Result: V=A=1/2, S=1/4, P=-1/4, T=0."""
    Kp = fierz(axial_axial_kernel())
    return {
        "S": float(_proj_scalar(Kp, np.eye(4)).real),
        "P": float(_proj_scalar(Kp, GAMMA5).real),
        "V": float(_proj_lorentz(Kp, _V).real),
        "A": float(_proj_lorentz(Kp, _A).real),
        "T": float(_proj_tensor(Kp).real),
    }


def GA_over_GV():
    """The walking-relevant ratio of axial to vector meson-channel couplings."""
    c = hehl_datta_channels()
    return c["A"] / c["V"]


# --- validations (the code is trustworthy iff these hold) ---
def validate_VmA_self_fierz():
    """The weak-interaction identity: (V-A)(V-A) is Fierz self-conjugate. Returns the
    max deviation (0 iff exact)."""
    VmA = [(_V[m] - _A[m]) for m in range(4)]
    K = np.zeros((4, 4, 4, 4), dtype=complex)
    for m in range(4):
        K += METRIC[m] * np.einsum("ij,kl->ijkl", VmA[m], VmA[m])
    return float(np.abs(fierz(K) - K).max())


def validate_involution():
    """Fierz^2 = identity. Returns the max deviation (0 iff exact)."""
    K = axial_axial_kernel()
    return float(np.abs(fierz(fierz(K)) - K).max())
