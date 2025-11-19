import numpy as np

def reaction_intensity(w0, sigma, h, Mf, fuel_params):
    """Reaction intensity (kW/m²) – Albini 1976."""
    if isinstance(w0, dict):
        p = w0
        w0, sigma, h = p["w0"], p["sigma"], p["h"]
    Gamma_prime = 0.0005 * sigma ** 1.5
    Gamma = Gamma_prime * ((192 + 0.2595 * sigma) ** -0.5)
    eta_M = 1 - 2.59*(Mf/0.3) + 5.11*(Mf/0.3)**2 - 3.52*(Mf/0.3)**3
    eta_M = np.clip(eta_M, 0.0, 1.0)
    Qig = 250 + 1116 * Mf
    return Gamma * w0 * h * eta_M / Qig * 1000  # kW/m²

def wind_factor(U, sigma):
    """Wind factor φ_w – Rothermel 1972."""
    C = 7.47 * np.exp(-0.133 * sigma ** 0.55)
    B = 0.02526 * sigma ** 0.54
    E = 0.715 * np.exp(-3.59e-4 * sigma)
    return C * (U ** B) * ((sigma / 3000) ** E)

def slope_factor(tan_phi):
    """Slope factor φ_s."""
    return 5.275 * (tan_phi ** 2)

def rothermel_rate_of_spread(fuel_model, Mf, U, tan_phi, fuel_params):
    """
    Forward rate of spread (m/s) – Rothermel 1972.

    Parameters
    ----------
    fuel_model : int
    Mf : float
        Dead fuel moisture (fraction)
    U : float
        Midflame wind speed (m/s)
    tan_phi : float
        Slope steepness (rise/run)
    fuel_params : dict

    Returns
    -------
    R : float
        Rate of spread (m/s)
    """
    p = fuel_params[fuel_model]
    w0 = p["w0"] * 20.83  # convert tons/acre to kg/m²
    sigma = p["sigma"]
    h = p["h"] * 2326  # BTU/lb → kJ/kg

    Ir = reaction_intensity(w0, sigma, h, Mf, None)
    phi_w = wind_factor(U, sigma)
    phi_s = slope_factor(tan_phi)
    xi = np.exp((0.792 + 0.681 * np.sqrt(sigma)) * (0.1 + phi_w + phi_s)) / (192 + 0.2595 * sigma)

    rho_b = w0 / p["delta"]  # bulk density kg/m³
    epsilon = 1 / (np.exp(0.43 * sigma) + 1)  # effective heating number
    Qig = 250 + 1116 * Mf

    R = (Ir * xi * (1 + phi_w + phi_s)) / (rho_b * epsilon * Qig)
    return R

def byram_flame_length(I_r):
    """Byram (1959) flame length (m)."""
    return 0.0775 * I_r ** 0.46
