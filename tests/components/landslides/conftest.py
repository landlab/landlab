import pytest

import numpy as np
from landlab import RasterModelGrid
from landlab.components import LandslideProbability


@pytest.fixture
def ls_prob():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__slope", at="node", dtype=float)
    grid.add_zeros("topographic__specific_contributing_area", at="node")
    grid.add_zeros("soil__transmissivity", at="node")
    grid.add_zeros("soil__saturated_hydraulic_conductivity", at="node")
    grid.add_zeros("soil__mode_total_cohesion", at="node")
    grid.add_zeros("soil__minimum_total_cohesion", at="node")
    grid.add_zeros("soil__maximum_total_cohesion", at="node")
    grid.add_zeros("soil__internal_friction_angle", at="node")
    grid.add_zeros("soil__density", at="node")
    grid.add_zeros("soil__thickness", at="node")

    return LandslideProbability(grid)


@pytest.fixture
def example_raster_model_grid():
    grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    gridnum = grid.number_of_nodes
    np.random.seed(seed=7)

    grid.at_node["topographic__slope"] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum).astype(float)
    grid.at_node["topographic__specific_contributing_area"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid.at_node["soil__saturated_hydraulic_conductivity"] = np.sort(
        np.random.randint(5, 20, gridnum).astype(float), -1
    )    
    hs = np.sort(np.round(100*np.random.uniform(1, 2, gridnum).astype(float))/100
    )
    grid.at_node["soil__thickness"] = hs
    grid.at_node["soil__transmissivity"] = \
    (
     grid.at_node["soil__saturated_hydraulic_conductivity"]*
     grid.at_node["soil__thickness"] 
    )    
    grid.at_node["soil__mode_total_cohesion"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid.at_node["soil__minimum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] - scatter_dat
    )
    grid.at_node["soil__maximum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] + scatter_dat
    )
    grid.at_node["soil__internal_friction_angle"] = np.sort(
        np.random.randint(26, 37, gridnum).astype(float)
    )

    grid.at_node["soil__density"] = 2000.0 * np.ones(gridnum)
    
    sat_thick = np.sort(np.round(100*np.random.uniform(1, 1.5, gridnum).astype(float))/100
    )        
    sat_thick[sat_thick>hs] = hs[sat_thick>hs]    
    
    grid.at_node["saturated__thickness"] = sat_thick    
    
    grid.at_node["mean__saturated_thickness"] = np.sort(np.round(100*np.random.uniform(0.75, 1.1, gridnum).astype(float))/100
    )
    
    grid.at_node["stdev__saturated_thickness"] = np.sort(np.round(100*np.random.uniform(0.1, 0.3, gridnum).astype(float))/100
    )
      
    return (grid)

