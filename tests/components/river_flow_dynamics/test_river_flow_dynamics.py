"""
Unit tests for landlab.components.river_flow_dynamics.river_flow_dynamics

last updated: 10/15/2023
"""
import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import river_flow_dynamics


@pytest.fixture
def rfd():
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    return river_flow_dynamics(grid)


def test_name(rfd):
    """Test component name"""
    print(rfd.name)
    assert rfd.name == "river_flow_dynamics"


def test_input_var_names(rfd):
    """Test input variable names"""
    assert rfd.input_var_names == (
        "surface_water__depth",
        "surface_water__elevation",
        "surface_water__velocity",
        "topographic__elevation",
    )


def test_output_var_names(rfd):
    """Test output variable names"""
    assert rfd.output_var_names == (
        "surface_water__depth",
        "surface_water__elevation",
        "surface_water__velocity",
    )


def test_var_units(rfd):
    assert rfd.var_units("surface_water__depth") == "m"
    assert rfd.var_units("surface_water__elevation") == "m"
    assert rfd.var_units("surface_water__velocity") == "m/s"
    assert rfd.var_units("topographic__elevation") == "m"


def test_initialization():
    """Test initialization with various parameters."""
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    rfd = river_flow_dynamics(grid)

    # Making sure fields have been created
    for field_name in rfd._info:
        if rfd._info[field_name]["mapping"] == "node":
            assert field_name in rfd.grid.at_node
        elif rfd._info[field_name]["mapping"] == "link":
            assert field_name in rfd.grid.at_link

    # Re-initialize, this time with fields already existing in the grid
    # (this triggers the "if" instead of "else" in the field setup in init)
    rfd = river_flow_dynamics(grid)


def test_run_one_step():
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    # We set fixed boundary conditions, specifying the nodes and links in which
    # the water is flowing into the grid
    fixed_entry_nodes = grid.nodes_at_left_edge
    fixed_entry_links = grid.links_at_node[fixed_entry_nodes][:, 0]

    # We set the fixed values in the entry nodes/links
    entry_nodes_h_values = 0.5 * np.ones_like(fixed_entry_nodes)
    entry_links_vel_values = 0.5 * np.ones_like(fixed_entry_links)

    # We specify the time step duration and we run it
    dt = 0.1
    rfd = river_flow_dynamics(
        grid,
        dt=dt,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    for _ in range(100):
        rfd.run_one_step()

    water_depth_solution = np.array(
        [
            0.4357753,
            0.4357753,
            0.43611027,
            0.43624251,
            0.43626605,
            0.43595278,
            0.43534349,
            0.43491662,
            0.43342158,
            0.43342158,
        ]
    )
    water_depth_obtained = grid.at_node["surface_water__depth"][
        grid.nodes_at_right_edge
    ]
    np.testing.assert_array_almost_equal(
        water_depth_solution, water_depth_obtained, decimal=6
    )
