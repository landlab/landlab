# landlab/components/fire_spread/utils.py

import numpy as np


def reaction_intensity(w0, sigma, h, Mf, _=None):
    """Reaction intensity Ir (kW/m²) – Albini version used in Rothermel."""
    # Moisture damping coefficient
    Mr = Mf / 0.30
    eta_M = 1.0 - 2.59 * Mr + 5.11 * Mr**2 - 3.52 * Mr**3
    eta_M = np.clip(eta_M, 0.0, 1.0)

    # Optimum reaction velocity
    Gamma_max = np.maximum(0.0592 * sigma**1.116, 1e-12)
    Gamma = Gamma_max * eta_M

    Qig = 250.0 + 1116.0 * Mf
    Ir = Gamma * w0 * h * eta_M / Qig * 1000.0  # kW/m²
    return Ir


def wind_factor(U, sigma):
    """Wind factor φ_w."""
    C = 7.47 * np.exp(-0.133 * sigma**0.55)
    B = 0.02526 * sigma**0.54
    E = 0.715 * np.exp(-0.000359 * sigma)
    phi_w = C * (U**B) * (sigma / 3000.0) ** (-E)
    return phi_w


def slope_factor(tan_phi):
    """Slope factor φ_s."""
    return 5.275 * tan_phi**2


def rothermel_rate_of_spread(fuel_model, Mf, U, tan_phi, fuel_params):
    """
    Rothermel (1972) forward rate of spread (m/s).
    Numerically stable implementation.
    """
    p = fuel_params[fuel_model]

    # Convert units
    w0 = p["w0"] * 20.83  # tons/acre → kg/m²
    sigma = p["sigma"]  # 1/ft  (already in correct units in Anderson 13)
    h = p["h"] * 2326.0  # BTU/lb → kJ/kg
    depth = p["delta"]  # fuel bed depth (m)

    # Reaction intensity
    Ir = reaction_intensity(w0, sigma, h, Mf)

    # Wind and slope factors
    phi_w = wind_factor(U, sigma)
    phi_s = slope_factor(tan_phi)

    # Propagating flux ratio ξ
    beta = w0 / (depth * 192.0 + 0.2595 * sigma)  # packing ratio
    beta_op = 3.348 * sigma ** (-0.8189)
    xi = np.exp((0.792 + 0.681 * np.sqrt(sigma)) * (beta_op + 0.1)) / (
        192.0 + 0.2595 * sigma
    )

    # Effective heating number ε
    epsilon = np.exp(-138.0 / sigma)

    # Bulk density
    rho_b = w0 / depth

    # Heat of pre-ignition
    Qig = 250.0 + 1116.0 * Mf

    # Final ROS
    denominator = rho_b * epsilon * Qig
    numerator = Ir * xi * (1.0 + phi_w + phi_s)

    R = numerator / np.where(denominator > 0, denominator, 1e-12)
    return np.clip(R, 0.0, None)  # never negative


def byram_flame_length(I_r):
    """Byram (1959) flame length (m)."""
    return 0.0775 * I_r**0.46
