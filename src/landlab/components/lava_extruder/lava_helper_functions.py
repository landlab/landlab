import numpy as np

def calc_lava_speed(
    lava_thickness, surface_slope, *, lava_density, gravity, viscosity, yield_stress
):
    """Calculate lava speed."""
    lava_thickness = np.asarray(lava_thickness)
    surface_slope = np.asarray(surface_slope)
    lava_speed = np.maximum(
        lava_density
        * gravity
        / (3.0 * viscosity)
        * lava_thickness**2
        * np.abs(surface_slope)
        - yield_stress / (2.0 * viscosity) * lava_thickness,
        0.0,
    )
    return lava_speed
