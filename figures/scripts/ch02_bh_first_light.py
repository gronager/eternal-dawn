#!/usr/bin/env python3
r"""Phase-0 "first light": the universe's measured numbers as black-hole properties.

Computes the consistency checks that anchor the "universe = black-hole interior"
thesis, prints a results table, writes machine-readable output to
sims/output/, and renders figures/pdf/bh_first_light.pdf.

Run directly (with the sims venv) or via `make figures`.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import blackhole as bh
from cartasis_sims import constants as k

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def compute() -> dict:
    M = bh.mass_within_hubble(Omega=1.0)
    R_H = bh.hubble_radius()
    res = {
        "H0_km_s_Mpc": k.H0_kms_Mpc,
        "R_H_Gly": R_H / k.Gly,
        "rho_crit_kg_m3": bh.critical_density(),
        "M_hubble_kg": M,
        "M_hubble_Msun": M / k.M_sun,
        "R_s_over_R_H": bh.rs_over_rh(Omega=1.0),
        "Omega_tot": k.Omega_tot,
        "a_star_max": bh.max_spin_parameter(Omega=1.0),
        "T_hawking_K": bh.hawking_temperature(M),
        "T_gibbons_hawking_K": bh.gibbons_hawking_temperature(),
        "cmb_over_T_hawking": k.T_cmb / bh.hawking_temperature(M),
        "filter_fraction_f": bh.filter_fraction(),
        "dm_to_baryon": k.Omega_c_hsq / k.Omega_b_hsq,
    }
    return res


def print_table(r: dict) -> str:
    lines = [
        "Phase-0: the universe's measured numbers as black-hole properties",
        "=" * 66,
        f"  Hubble radius           R_H   = {r['R_H_Gly']:.2f} Gly",
        f"  Critical density        rho_c = {r['rho_crit_kg_m3']:.3e} kg/m^3",
        f"  Hubble-sphere mass       M     = {r['M_hubble_kg']:.3e} kg "
        f"({r['M_hubble_Msun']:.2e} Msun)",
        "",
        f"  [A] R_s / R_H = Omega          = {r['R_s_over_R_H']:.6f}   "
        "(flat universe sits at its Schwarzschild radius)",
        f"  [B] max causal spin a*_max     = {r['a_star_max']:.4f}     "
        "(order unity -> realistic spin is sub-extremal)",
        f"  [C] Hawking T of M             = {r['T_hawking_K']:.3e} K",
        f"      Gibbons-Hawking T (now)    = {r['T_gibbons_hawking_K']:.3e} K "
        "(= 2 x Hawking T)",
        f"      CMB / Hawking-T            = {r['cmb_over_T_hawking']:.2e}   "
        "(30 orders: CMB is NOT present-horizon radiation)",
        f"  [D] filter pass-fraction f     = {r['filter_fraction_f']:.4f}    "
        f"(DM/baryon = {r['dm_to_baryon']:.2f} -> f ~ 1/{1/r['filter_fraction_f']:.1f})",
        "=" * 66,
    ]
    return "\n".join(lines)


def make_figure(r: dict) -> str:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel 1: Hawking temperature vs mass, with landmarks.
    M = np.logspace(0, 54, 400)  # kg
    T = bh.hawking_temperature(M)
    ax1.loglog(M, T, color="0.2", lw=1.6)
    ax1.axhline(k.T_cmb, color="C3", ls="--", lw=1.2)
    ax1.text(2e1, k.T_cmb * 1.6, f"CMB  {k.T_cmb} K", color="C3", fontsize=9)

    landmarks = {
        "stellar BH\n(10 M$_\\odot$)": 10 * k.M_sun,
        "SMBH\n(10$^9$ M$_\\odot$)": 1e9 * k.M_sun,
        "observable\nuniverse": r["M_hubble_kg"],
    }
    for label, mass in landmarks.items():
        ax1.plot([mass], [bh.hawking_temperature(mass)], "o", color="C0", ms=6)
        ax1.annotate(label, (mass, bh.hawking_temperature(mass)),
                     textcoords="offset points", xytext=(-4, 8),
                     fontsize=8, ha="center")
    ax1.set_xlabel("black-hole mass  $M$  [kg]")
    ax1.set_ylabel("Hawking temperature  $T_H$  [K]")
    ax1.set_title("Hawking temperature vs mass")
    ax1.grid(True, which="major", alpha=0.25)

    # Panel 2: the dimensionless identities as a bar chart vs their "1" scale.
    names = [r"$R_s/R_H=\Omega$", r"$a^*_{\max}$", r"filter $f$"]
    vals = [r["R_s_over_R_H"], r["a_star_max"], r["filter_fraction_f"]]
    refs = [1.0, 1.0, None]
    colors = ["C2", "C0", "C1"]
    xpos = np.arange(len(names))
    ax2.bar(xpos, vals, color=colors, width=0.55, alpha=0.85)
    ax2.axhline(1.0, color="0.5", ls=":", lw=1.0)
    ax2.text(len(names) - 0.5, 1.03, "extremal / unity", color="0.4",
             fontsize=8, ha="right")
    for x, v in zip(xpos, vals):
        ax2.text(x, v + 0.05, f"{v:.3g}", ha="center", fontsize=9)
    ax2.set_xticks(xpos)
    ax2.set_xticklabels(names)
    ax2.set_ylim(0, 2.3)
    ax2.set_title("Black-hole identities from measured cosmology")

    fig.suptitle("Universe as a black-hole interior: Phase-0 consistency checks",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(PDF_DIR, "bh_first_light.pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    return out


def main() -> None:
    r = compute()
    table = print_table(r)
    print(table)
    with open(os.path.join(OUT_DIR, "first_light.txt"), "w") as f:
        f.write(table + "\n")
    with open(os.path.join(OUT_DIR, "first_light.json"), "w") as f:
        json.dump(r, f, indent=2)
    pdf = make_figure(r)
    print(f"\nwrote {pdf}")
    print(f"wrote {os.path.join(OUT_DIR, 'first_light.txt')}")


if __name__ == "__main__":
    main()
