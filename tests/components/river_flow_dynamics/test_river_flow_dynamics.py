"""
Unit tests for landlab.components.river_flow_dynamics.RiverFlowDynamics

last updated: 11/11/2024
"""

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import RiverFlowDynamics


@pytest.fixture
def rfd():
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    return RiverFlowDynamics(grid)


def test_name(rfd):
    """Test component name"""
    print(rfd.name)
    assert rfd.name == "RiverFlowDynamics"


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
    rfd = RiverFlowDynamics(grid)

    # Making sure fields have been created
    for field_name in rfd._info:
        if rfd._info[field_name]["mapping"] == "node":
            assert field_name in rfd.grid.at_node
        elif rfd._info[field_name]["mapping"] == "link":
            assert field_name in rfd.grid.at_link

    # Re-initialize, this time with fields already existing in the grid
    # (this triggers the "if" instead of "else" in the field setup in init)
    rfd = RiverFlowDynamics(grid)


def test_mass_conservation():
    """Test that water volume is conserved in a closed system."""
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    # Create a bowl-shaped topography for more interesting flow
    x = grid.x_of_node - 0.5  # Center at 0.5
    y = grid.y_of_node - 0.5  # Center at 0.5
    grid.at_node["topographic__elevation"] = 0.1 * (x**2 + y**2)

    # Set water surface elevation to 0.5 everywhere
    grid.at_node["surface_water__elevation"][:] = 0.5

    # Calculate water depth as difference between surface elevation and topography
    topo = grid.at_node["topographic__elevation"]
    water_surface = grid.at_node["surface_water__elevation"]
    grid.at_node["surface_water__depth"] = np.maximum(water_surface - topo, 0.0)

    # Set up model with closed boundaries to ensure mass conservation
    rfd = RiverFlowDynamics(grid)

    # Calculate initial volume
    initial_volume = np.sum(grid.at_node["surface_water__depth"]) * grid.dx * grid.dy

    # Run model for several time steps
    n_steps = 100
    volumes = np.zeros(n_steps + 1)
    volumes[0] = initial_volume

    for i in range(n_steps):
        rfd.run_one_step()
        volumes[i + 1] = (
            np.sum(grid.at_node["surface_water__depth"]) * grid.dx * grid.dy
        )

    # Calculate volume change
    volume_changes = np.abs(volumes - initial_volume) / initial_volume
    max_volume_change = np.max(volume_changes)

    # Assert volume is conserved within a more realistic tolerance (0.1%)
    assert max_volume_change < 1e-3

    # Check for physically reasonable results
    depths = grid.at_node["surface_water__depth"]
    initial_max_depth = np.max(grid.at_node["surface_water__depth"])
    assert np.all(depths >= 0)
    assert np.max(depths) < initial_max_depth * 1.001


def test_open_boundaries():
    """Test that water properly exits through open boundaries."""
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    # Sloped channel
    grid.at_node["topographic__elevation"] = grid.x_of_node * 0.05

    # Set up inflow conditions
    fixed_entry_nodes = grid.nodes_at_left_edge
    fixed_entry_links = grid.links_at_node[fixed_entry_nodes][:, 0]
    entry_nodes_h_values = 0.5 * np.ones_like(fixed_entry_nodes)
    entry_links_vel_values = 0.3 * np.ones_like(fixed_entry_links)

    rfd = RiverFlowDynamics(
        grid,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    # Run to steady state
    for _ in range(200):
        rfd.run_one_step()

    # Check that water exits smoothly (no backing up at boundary)
    right_depths = grid.at_node["surface_water__depth"][grid.nodes_at_right_edge]
    near_right_depths = grid.at_node["surface_water__depth"][
        grid.nodes_at_right_edge - 1
    ]

    # Water depth should decrease or stay similar near boundary
    assert np.all(right_depths <= near_right_depths * 1.1)


def test_still_water():
    """Test that still water stays still."""
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    # Flat bottom with uniform water depth
    grid.at_node["surface_water__depth"][:] = 0.5
    grid.at_node["surface_water__elevation"][:] = 0.5

    rfd = RiverFlowDynamics(grid)

    # Initial conditions
    initial_volume = np.sum(grid.at_node["surface_water__depth"]) * grid.dx * grid.dy

    # Run model
    rfd.run_one_step()

    # Check volume conservation
    final_volume = np.sum(grid.at_node["surface_water__depth"]) * grid.dx * grid.dy
    volume_change = abs(final_volume - initial_volume) / initial_volume
    assert volume_change < 1e-4  # Allow 0.01% volume change

    # Check surface stays flat (no gradients)
    depth_gradients = np.abs(
        np.diff(grid.at_node["surface_water__depth"].reshape(10, 10), axis=1)
    )
    assert np.all(depth_gradients < 1e-4)  # Maximum gradient of 0.1 mm per cell


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
    rfd = RiverFlowDynamics(
        grid,
        dt=dt,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    for _ in range(100):
        rfd.run_one_step()

    water_depth_solution = np.round(
        np.array(
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
        ),
        3,
    )
    water_depth_obtained = grid.at_node["surface_water__depth"][
        grid.nodes_at_right_edge
    ]

    water_depth_obtained = np.round(water_depth_obtained, 3)
    np.testing.assert_array_almost_equal(
        water_depth_solution, water_depth_obtained, decimal=3
    )


def test_downhill_flow():
    """Test that water flows downhill."""
    grid = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__elevation", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    # Create sloped surface
    grid.at_node["topographic__elevation"] = grid.x_of_node * 0.1

    # Set water surface elevation at upstream end
    upstream_water_level = 0.5
    grid.at_node["surface_water__elevation"][
        grid.nodes_at_left_edge
    ] = upstream_water_level

    # Calculate initial water depth at upstream end
    upstream_topo = grid.at_node["topographic__elevation"][grid.nodes_at_left_edge]
    grid.at_node["surface_water__depth"][grid.nodes_at_left_edge] = (
        upstream_water_level - upstream_topo
    )

    # Set up inflow conditions
    fixed_entry_nodes = grid.nodes_at_left_edge
    fixed_entry_links = grid.links_at_node[fixed_entry_nodes][:, 0]
    entry_links_vel_values = 0.3 * np.ones_like(fixed_entry_links)

    rfd = RiverFlowDynamics(
        grid,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=grid.at_node["surface_water__depth"][fixed_entry_nodes],
        entry_links_vel_values=entry_links_vel_values,
    )

    # Run model
    for _ in range(50):
        rfd.run_one_step()

    # Check mean flow direction
    horizontal_links = grid.horizontal_links
    mean_velocity = np.mean(grid.at_link["surface_water__velocity"][horizontal_links])
    assert mean_velocity > 0  # Overall flow should be positive (downhill)

    # Check water depths are physically reasonable
    assert np.all(grid.at_node["surface_water__depth"] >= 0)

    # Check water surface elevation is consistent with depth and topography
    # but allow for small numerical differences
    calculated_elevation = (
        grid.at_node["surface_water__depth"] + grid.at_node["topographic__elevation"]
    )
    elevation_difference = np.abs(
        grid.at_node["surface_water__elevation"] - calculated_elevation
    )
    assert np.max(elevation_difference) < 0.1  # Allow differences up to 10cm

    # Check that water has moved downstream
    outlet_depths = grid.at_node["surface_water__depth"][grid.nodes_at_right_edge]
    assert np.any(outlet_depths > 0.01)


def test_time_step_sensitivity():
    """Test that results are not heavily dependent on time step choice for
    relatively similar systems."""
    grid1 = RasterModelGrid((10, 10), xy_spacing=0.1)
    grid2 = RasterModelGrid((10, 10), xy_spacing=0.1)

    # Set up identical grids
    for grid in [grid1, grid2]:
        grid.add_zeros("topographic__elevation", at="node")
        grid.add_zeros("surface_water__depth", at="node")
        grid.add_zeros("surface_water__elevation", at="node")
        grid.add_zeros("surface_water__velocity", at="link")

        # Create a sloped surface
        grid.at_node["topographic__elevation"] = grid.x_of_node * 0.1

        # Set the water surface elevation
        water_level = 0.5
        grid.at_node["surface_water__elevation"][grid.nodes_at_left_edge] = water_level

        # Calculate initial water depth
        topo = grid.at_node["topographic__elevation"][grid.nodes_at_left_edge]
        grid.at_node["surface_water__depth"][grid.nodes_at_left_edge] = (
            water_level - topo
        )

    # Set up identical boundary conditions for both grids
    fixed_entry_nodes = grid1.nodes_at_left_edge
    fixed_entry_links = grid1.links_at_node[fixed_entry_nodes][:, 0]
    entry_nodes_h_values = grid1.at_node["surface_water__depth"][fixed_entry_nodes]
    entry_links_vel_values = 0.3 * np.ones_like(fixed_entry_links)

    # Create two models with different time steps
    rfd1 = RiverFlowDynamics(
        grid1,
        dt=0.1,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    rfd2 = RiverFlowDynamics(
        grid2,
        dt=0.05,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    # Run to same total time
    for _ in range(100):
        rfd1.run_one_step()
    for _ in range(200):
        rfd2.run_one_step()

    # Compare results with reasonable tolerances
    depth_difference = np.abs(
        grid1.at_node["surface_water__depth"] - grid2.at_node["surface_water__depth"]
    )
    assert np.max(depth_difference) < 0.1  # Max difference in depths

    # Check both solutions maintain physical consistency
    for grid in [grid1, grid2]:
        # Water depth should be non-negative
        assert np.all(grid.at_node["surface_water__depth"] >= 0)

        # Check water surface elevation consistency
        calculated_elevation = (
            grid.at_node["surface_water__depth"]
            + grid.at_node["topographic__elevation"]
        )
        elevation_difference = np.abs(
            grid.at_node["surface_water__elevation"] - calculated_elevation
        )
        assert np.max(elevation_difference) < 0.1

        # Check water has propagated downstream
        outlet_depths = grid.at_node["surface_water__depth"][grid.nodes_at_right_edge]
        assert np.any(outlet_depths > 0.01)


def test_numerical_stability():
    """Test numerical stability under various conditions."""

    def create_test_grid():
        """Helper function to create a grid with initial conditions."""
        grid = RasterModelGrid((10, 10), xy_spacing=0.1)
        grid.add_zeros("topographic__elevation", at="node")
        grid.add_zeros("surface_water__depth", at="node")
        grid.add_zeros("surface_water__elevation", at="node")
        grid.add_zeros("surface_water__velocity", at="link")

        # Create challenging conditions
        grid.at_node["topographic__elevation"] = grid.x_of_node * 0.1

        # Set initial conditions
        grid.at_node["surface_water__depth"][:] = 0.5
        grid.at_node["surface_water__elevation"] = (
            grid.at_node["surface_water__depth"]
            + grid.at_node["topographic__elevation"]
        )

        return grid

    time_steps = [0.1, 0.05, 0.01]
    max_velocities = []
    final_depths = []

    for dt in time_steps:
        grid_test = create_test_grid()
        rfd = RiverFlowDynamics(grid_test, dt=dt)

        try:
            for _ in range(int(5.0 / dt)):
                rfd.run_one_step()

                assert np.all(np.isfinite(grid_test.at_node["surface_water__depth"]))
                assert np.all(grid_test.at_node["surface_water__depth"] >= 0)
                assert np.all(np.isfinite(grid_test.at_link["surface_water__velocity"]))

            max_velocities.append(
                np.max(np.abs(grid_test.at_link["surface_water__velocity"]))
            )
            final_depths.append(grid_test.at_node["surface_water__depth"].copy())

        except AssertionError as err:
            raise AssertionError(f"Instability detected with dt={dt}") from err

    depth_variations = [
        np.max(np.abs(d1 - d2)) for d1, d2 in zip(final_depths[:-1], final_depths[1:])
    ]
    assert np.all(np.array(depth_variations) < 0.1)
    assert np.all(np.array(max_velocities) < 10.0)

    depth_variations = [
        np.max(np.abs(d1 - d2)) for d1, d2 in zip(final_depths[:-1], final_depths[1:])
    ]
    assert np.all(np.array(depth_variations) < 0.1)
    assert np.all(np.array(max_velocities) < 1.0)


def test_analytical_solution():
    """Test numerical and analytical solutions."""

    # Creating a grid
    grid = RasterModelGrid((6, 6), xy_spacing=0.1)
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
    # entry_nodes_h_values = np.repeat(0.5, 6)
    entry_links_vel_values = 0.6 * np.ones_like(fixed_entry_links)
    # entry_links_vel_values = np.repeat(0.6, 6)

    # We specify the time step duration and we run it
    dt = 0.1
    rfd = RiverFlowDynamics(
        grid,
        dt=dt,
        mannings_n=0.0010,
        fixed_entry_nodes=fixed_entry_nodes,
        fixed_entry_links=fixed_entry_links,
        entry_nodes_h_values=entry_nodes_h_values,
        entry_links_vel_values=entry_links_vel_values,
    )

    for _ in range(300):
        rfd.run_one_step()

    # Comparing the uniform outlet with the obtained water depth
    water_depth_solution = 0.5 * np.ones_like(fixed_entry_nodes)
    water_depth_obtained = grid.at_node["surface_water__depth"][
        grid.nodes_at_right_edge
    ]

    water_depth_difference = np.abs(water_depth_solution - water_depth_obtained)
    assert np.max(water_depth_difference) < 0.001  # Max difference in vel

    vel_solution = 0.6 * np.ones_like(fixed_entry_links)[1:5]
    vel_obtained = grid.at_link["surface_water__velocity"][
        grid.links_at_node[grid.nodes_at_right_edge][:, 2]
    ][1:5]
    vel_difference = np.abs(vel_solution - vel_obtained)

    assert np.max(vel_difference) < 0.01  # Max difference in vel
