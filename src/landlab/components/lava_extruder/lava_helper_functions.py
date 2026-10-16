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


def conductive_heat_flux(
    lava_temperature,
    area,
    conductive_thickness,
    *,
    thermal_conductivity,
    ambient_temperature,
):
    """Calculate conductive heat flux."""
    lava_temperature = np.asarray(lava_temperature)
    area = np.asarray(area)
    conductive_thickness = np.asarray(conductive_thickness)

    heat_flux = (
        thermal_conductivity
        * (lava_temperature - ambient_temperature)
        / conductive_thickness
        * area
    )
    return heat_flux


def thermal_energy_to_temperature(
    energy, lava_thickness, cell_area, *, lava_density, specific_heat_capacity
):
    """Convert thermal energy to temperature."""
    energy = np.asarray(energy)
    lava_thickness = np.asarray(lava_thickness)
    cell_area = np.asarray(cell_area)

    temperature = energy / (
        lava_thickness * cell_area * lava_density * specific_heat_capacity
    )
    return temperature


def vent_thermal_energy(
    temperature, lava_thickness, cell_area, *, lava_density, specific_heat_capacity
):
    """Convert temperature to thermal energy."""
    temperature = np.asarray(temperature)
    lava_thickness = np.asarray(lava_thickness)
    cell_area = np.asarray(cell_area)

    energy = (
        temperature * lava_thickness * cell_area * lava_density * specific_heat_capacity
    )
    return energy
