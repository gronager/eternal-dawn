#!/usr/bin/env python3
r"""The tensor-to-scalar ratio r: the inflation-vs-bounce fork (Ch.7).

Left -- the observational (n_s, r) landscape: the Planck n_s band, the current
BICEP/Keck upper bound r < 0.036, LiteBIRD's forecast reach (sigma_r ~ 0.001), and
the inflation menu (Starobinsky/Higgs plateau at r ~ 0.0037; m^2 phi^2 large-field
at r ~ 0.14, already excluded). The matter/torsion bounce sits at n_s ~ 0.965 with r
NOT fixed by a consistency relation: the simplest single-field version over-produces
tensors (r >~ 1, excluded) and must be cured by enhancing the scalar power, which
generically drops r below LiteBIRD.

Right -- the discriminator that does not depend on r's uncertain magnitude: the
tensor TILT. Inflation ties n_T to r (n_T = -r/8, so plateau models sit at n_T ~ 0);
the matter bounce ties it to the SCALAR tilt (n_T = n_s - 1 ~ -0.035). At the
observed n_s the two predictions for n_T differ ~70x -- a clean, in-principle
measurable fork.

Renders figures/pdf/tensor_r.pdf, writes sims/output/tensor_r.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import tensor_ratio as tr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    ns_star, r_star = tr.starobinsky(57.0)
    ns_phi2, r_phi2 = tr.chaotic_phi2(57.0)
    nT_bounce = tr.bounce_tensor_tilt(tr.N_S_OBS)
    nT_infl = tr.inflation_consistency_nT(r_star)
    d = tr.discriminator_gap()
    tilt_ratio = abs(nT_bounce / nT_infl)   # how many x more tilt the bounce predicts

    lines = [
        "The tensor-to-scalar ratio r: inflation vs the matter/torsion bounce",
        "=" * 68,
        f"  Planck n_s            = {tr.N_S_OBS:.4f} +/- {tr.N_S_SIGMA:.4f}",
        f"  current bound         = r < {tr.BK18_UPPER:.3f}  (BICEP/Keck 2021)",
        f"  LiteBIRD reach        = sigma(r) ~ {tr.LITEBIRD_SIGMA_R:.0e}",
        "",
        "  inflation menu (N=57 e-folds):",
        f"    Starobinsky/Higgs   (n_s, r) = ({ns_star:.4f}, {r_star:.4f})  "
        "<- within LiteBIRD reach",
        f"    m^2 phi^2 (chaotic) (n_s, r) = ({ns_phi2:.4f}, {r_phi2:.3f})  "
        "<- excluded (r > 0.036)",
        "",
        "  the bounce's robust predictions (primordial.py):",
        f"    n_s   = 0.965  (contraction just softer than matter, w ~ "
        f"{tr.w_for_observed_tilt():+.4f})",
        f"    n_T   = n_s - 1 = {nT_bounce:+.4f}   (equal-tilt: tensor tilt = scalar"
        " tilt)",
        "    r     : NOT fixed -- simplest single-field bounce gives r ~ O(1)"
        " (excluded),",
        "            cured (small c_s / 2nd field) -> generically below LiteBIRD.",
        "",
        "  the tilt discriminator at the observed n_s:",
        f"    bounce    n_T = {d['nT_bounce']:+.4f}",
        f"    inflation n_T = {d['nT_inflation']:+.4f}   (= -r/8, Starobinsky r)",
        f"    gap       = {d['gap']:+.4f}   (bounce tilt ~{tilt_ratio:.0f}x"
        " inflation's -- the clean fork)",
        "",
        "READING: r's magnitude alone is not yet a sharp bounce prediction (the",
        "simplest matter bounce over-produces tensors and must be cured). The robust",
        "fork is the CONSISTENCY RELATION: inflation puts (n_T, r) on n_T = -r/8;",
        "the bounce puts n_T = n_s - 1, decoupled from r. A plateau-value DETECTION on",
        "the inflationary line favours inflation; a NULL (r < ~0.001) or a measured",
        "n_T ~ n_s-1 favours the bounce. LiteBIRD decides.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "tensor_r.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.4, 5.4))

    # ---- left: (n_s, r) landscape -------------------------------------
    axL.axvspan(tr.N_S_OBS - tr.N_S_SIGMA, tr.N_S_OBS + tr.N_S_SIGMA,
                color="C2", alpha=0.15, label=r"Planck $n_s$")
    axL.axhspan(tr.BK18_UPPER, 0.25, color="0.85", alpha=0.7)
    axL.axhline(tr.BK18_UPPER, color="0.4", lw=1.2, ls="--")
    axL.text(0.943, tr.BK18_UPPER * 1.15, "excluded  ($r<0.036$, BK18)",
             fontsize=8.5, color="0.35")
    axL.axhline(tr.LITEBIRD_SIGMA_R, color="C4", lw=1.2, ls=":")
    axL.text(0.943, tr.LITEBIRD_SIGMA_R * 1.15, r"LiteBIRD reach $\sigma_r\sim10^{-3}$",
             fontsize=8.5, color="C4")
    # inflation models
    axL.plot(ns_star, r_star, "C0o", ms=9)
    axL.annotate("Starobinsky / Higgs\n(plateau)", (ns_star, r_star),
                 textcoords="offset points", xytext=(8, -2), fontsize=8.5, color="C0")
    axL.plot(ns_phi2, r_phi2, "rx", ms=10, mew=2)
    axL.annotate(r"$m^2\phi^2$ (excluded)", (ns_phi2, r_phi2),
                 textcoords="offset points", xytext=(-4, 8), fontsize=8.5, color="r")
    # the bounce band: n_s ~ 0.965, r a free/low output
    axL.axvline(tr.N_S_OBS, color="C3", lw=2.0)
    axL.annotate("matter / torsion bounce\n$n_s\\simeq0.965$, $r$ set by the bounce\n"
                 "(simplest $r\\gtrsim1$ excluded $\\to$ cured low)",
                 (tr.N_S_OBS, 3e-4), textcoords="offset points", xytext=(-150, 0),
                 fontsize=8.5, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=1.2))
    axL.annotate("", (tr.N_S_OBS, 2e-4), (tr.N_S_OBS, 8e-4),
                 arrowprops=dict(arrowstyle="<-", color="C3", lw=1.6))
    axL.set_yscale("log")
    axL.set_xlim(0.94, 0.99)
    axL.set_ylim(1e-4, 0.25)
    axL.set_xlabel(r"scalar tilt  $n_s$")
    axL.set_ylabel(r"tensor-to-scalar ratio  $r$")
    axL.set_title("The observational landscape", fontsize=11)
    axL.legend(fontsize=8.5, loc="lower left")
    axL.grid(True, which="both", alpha=0.2)

    # ---- right: the tilt discriminator (n_s, n_T) ---------------------
    ns = np.linspace(0.94, 0.99, 50)
    axR.plot(ns, ns - 1.0, "C3", lw=2.0, label=r"bounce: $n_T=n_s-1$")
    axR.axhline(nT_infl, color="C0", lw=2.0,
                label=r"inflation: $n_T=-r/8\approx0$ (plateau)")
    axR.axvspan(tr.N_S_OBS - tr.N_S_SIGMA, tr.N_S_OBS + tr.N_S_SIGMA,
                color="C2", alpha=0.15, label=r"Planck $n_s$")
    axR.plot(tr.N_S_OBS, nT_bounce, "C3o", ms=8)
    axR.plot(tr.N_S_OBS, nT_infl, "C0o", ms=8)
    axR.annotate("", (tr.N_S_OBS + 0.004, nT_bounce),
                 (tr.N_S_OBS + 0.004, nT_infl),
                 arrowprops=dict(arrowstyle="<->", color="0.3", lw=1.4))
    axR.text(tr.N_S_OBS + 0.006, 0.5 * (nT_bounce + nT_infl),
             f"$\\sim${tilt_ratio:.0f}$\\times$ gap\n"
             "(the fork)", fontsize=8.5, color="0.3", va="center")
    axR.set_xlim(0.94, 0.99)
    axR.set_ylim(-0.06, 0.02)
    axR.set_xlabel(r"scalar tilt  $n_s$")
    axR.set_ylabel(r"tensor tilt  $n_T$")
    axR.set_title("The consistency-relation fork", fontsize=11)
    axR.legend(fontsize=8.5, loc="lower right")
    axR.grid(True, alpha=0.2)

    fig.suptitle("Primordial gravitational waves: inflation vs the torsion bounce",
                 fontsize=12.5, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(PDF_DIR, "tensor_r.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "tensor_r.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'tensor_r.pdf')}")


if __name__ == "__main__":
    main()
