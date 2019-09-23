import pytest

from landlab import RasterModelGrid
from landlab.components.soil_moisture import SoilInfiltrationGreenAmpt
from landlab.components.soil_moisture.soil_moisture_dynamics import SoilMoisture


@pytest.fixture
def sm():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    return SoilMoisture(grid)


@pytest.fixture
def si():
    grid = RasterModelGrid((10, 10), xy_spacing=25)
    grid.add_ones("soil_water_infiltration__depth", at="node", dtype=float)
    grid.add_ones("surface_water__depth", at="node")
    hydraulic_conductivity = 2.5 * (10 ** -5)
    grid.at_node["surface_water__depth"] *= 0.5
    grid.at_node["soil_water_infiltration__depth"] *= 10 ** -5
    return SoilInfiltrationGreenAmpt(
        grid,
        hydraulic_conductivity=hydraulic_conductivity,
        soil_bulk_density=1700.0,
        rock_density=2650.0,
        initial_soil_moisture_content=0.2,
        soil_type="silt loam",
        volume_fraction_coarse_fragments=0.6,
        coarse_sed_flag=False,
        surface_water_minimum_depth=1.0e-7,
        soil_pore_size_distribution_index=None,
        soil_bubbling_pressure=None,
        wetting_front_capillary_pressure_head=None,
    )
