"""Physical constants (SI, CODATA 2018) and Planck 2018 cosmological parameters.

Everything in SI. Cosmology values are Planck 2018 TT,TE,EE+lowE+lensing
(Planck Collaboration VI 2020).
"""

from __future__ import annotations

# --- Fundamental constants (SI) ---
G = 6.67430e-11          # gravitational constant, m^3 kg^-1 s^-2
c = 2.99792458e8         # speed of light, m s^-1
hbar = 1.054571817e-34   # reduced Planck constant, J s
kB = 1.380649e-23        # Boltzmann constant, J K^-1
h_planck = 6.62607015e-34

# --- Astronomical units ---
Mpc = 3.0856775814913673e22   # metres
Gly = 9.4607304725808e24      # metres (10^9 light-years)
M_sun = 1.98892e30            # kg
year = 3.15576e7              # s

# --- Planck 2018 cosmological parameters ---
H0_kms_Mpc = 67.36       # km s^-1 Mpc^-1
h = 0.6736               # dimensionless Hubble parameter
Omega_m = 0.3153
Omega_Lambda = 0.6847
Omega_k = 0.0007         # |Omega_k| < 0.005: spatially flat to high precision
Omega_b_hsq = 0.02237    # physical baryon density  omega_b = Omega_b h^2
Omega_c_hsq = 0.1200     # physical cold dark matter density omega_c = Omega_c h^2
T_cmb = 2.7255           # K

# Derived
H0 = H0_kms_Mpc * 1.0e3 / Mpc          # s^-1
Omega_b = Omega_b_hsq / h**2
Omega_c = Omega_c_hsq / h**2
Omega_tot = Omega_m + Omega_Lambda     # ~1 by construction (flat)
