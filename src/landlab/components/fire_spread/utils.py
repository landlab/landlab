# landlab/components/fire_spread/utils.py

import numpy as np


def reaction_intensity(w0, sigma, h, Mf, _=None):
    """Reaction intensity Ir (kW/m²) – standard Albini/Rothermel form."""
    Mr = Mf / 0.30
    eta_M = 1.0 - 2.59*Mr + 5.11*Mr**2 - 3.52*Mr**3
    eta_M = np.clip(eta_M, 0.0, 1.0)

    Gamma_max = np.maximum(sigma**1.116 * 0.0592, 1e-12)
    Gamma = Gamma_max * eta_M

    Qig = 250.0 + 1116.0 * Mf
    Ir = Gamma * w0 * h * eta_M / Qig * 1000.0          # kW/m²
    return Ir


def wind_factor(U, sigma):
    C = 7.47 * np.exp(-0.133 * sigma**0.55)
    B = 0.02526 * sigma**0.54
    E  = 0.715 * np.exp(-0.000359 * sigma)
    return C * (U**B) * (sigma/3000.0)**(-E)


def slope_factor(tan_phi):
    return 5.275 * tan_phi**2


def rothermel_rate_of_spread(fuel_model, Mf, U, tan_phi, fuel_params):
    """
    Rothermel (1972) rate of spread in **meters per second**.
    Matches BehavePlus, FARSITE, and the official USFS implementation.
    """
    p = fuel_params[fuel_model]

    # Unit conversions
    w0    = p["w0"]   * 20.83      # tons/acre → kg/m²
    sigma = p["sigma"]             # already 1/ft in Anderson 13
    h     = p["h"]    * 2326.0     # BTU/lb → kJ/kg
    depth = p["delta"]             # feet → already in metres in Anderson 13

    Ir    = reaction_intensity(w0, sigma, h, Mf)
    phi_w = wind_factor(U, sigma)
    phi_s = slope_factor(tan_phi)

    # Propagating flux ratio
    beta_op = 3.348 * sigma**(-0.8189)
    xi = np.exp((0.792 + 0.681 * np.sqrt(sigma)) * (beta_op + 0.1)) / (192 + 0.2595*sigma)

    # Effective heating number
    epsilon = np.exp(-138.0 / sigma)

    # Packing ratio and bulk density
    rho_b = w0 / depth

    Qig = 250.0 + 1116.0 * Mf

    # Rate of spread in **feet per minute** (original Rothermel units)
    R_fpm = Ir * xi * (1 + phi_w + phi_s) / (rho_b * epsilon * Qig)

    # Convert to metres per second
    R = R_fpm * 0.00508   # 1 ft/min = 0.00508 m/s

    return np.clip(R, 0.0, None)


def byram_flame_length(I_r):
    """Byram flame length in metres."""
    return 0.0775 * I_r**0.46
