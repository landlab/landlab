"""
Unit tests for landlab.components.river_bed_dynamics.river_bed_dynamics

last updated: 10/12/2023
"""
import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import OverlandFlow, river_bed_dynamics
from landlab.grid.mappers import map_mean_of_link_nodes_to_link

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(r_b_d):
    assert r_b_d.name == "river_bed_dynamics"


def test_input_var_names(r_b_d):
    assert r_b_d.input_var_names == (
        "surface_water__depth",
        "surface_water__velocity",
        "topographic__elevation",
    )


def test_output_var_names(r_b_d):
    assert r_b_d.output_var_names == ("topographic__elevation",)


def test_optional_var_names(r_b_d):
    assert r_b_d.optional_var_names == ()


def test_var_units(r_b_d):
    assert set(r_b_d.input_var_names) | set(r_b_d.output_var_names) | set(
        r_b_d.optional_var_names
    ) == set(dict(r_b_d.units).keys())

    assert r_b_d.var_units("surface_water__depth") == "m"
    assert r_b_d.var_units("surface_water__velocity") == "m/s"
    assert r_b_d.var_units("topographic__elevation") == "m"


def test_grid_shape(r_b_d):
    assert r_b_d.grid.number_of_node_rows == _SHAPE[0]
    assert r_b_d.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(r_b_d):
    assert r_b_d.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(r_b_d):
    assert r_b_d.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_median_size():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [4, 10], [2, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    rbd = river_bed_dynamics(grid, gsd=gsd)
    assert_almost_equal(rbd._bed_surface__median_size_node[20], 16)


def test_geometric_mean_size_node():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [4, 10], [2, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    rbd = river_bed_dynamics(grid, gsd=gsd)
    np.testing.assert_almost_equal(
        rbd._bed_surface__geometric_standard_deviation_size_node[20], 2.567, decimal=3
    )


def test_sand_fraction_node():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [2, 10], [1, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    rbd = river_bed_dynamics(grid, gsd=gsd)
    np.testing.assert_almost_equal(
        rbd._bed_surface__sand_fraction_node[20], 0.1, decimal=3
    )


def test_velocity_previous_time():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    grid["link"]["surface_water__velocity"][20] = 10
    rbd = river_bed_dynamics(grid, gsd=gsd)

    np.testing.assert_almost_equal(
        rbd._surface_water__velocity_previous_time_link[20], 10, decimal=3
    )


def test_error_gsd_location_node():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="bed_surface__grain_size_distribution_location_node"
        ".*does not have the same dimensions of the grid's nodes",
    ):
        river_bed_dynamics(
            grid, gsd=gsd, bed_surface__grain_size_distribution_location_node=testArray
        )


def test_error_supply_imposed():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="sediment_transport__sediment_supply_imposed_link"
        ".*does not have the same dimensions of the grid's link",
    ):
        river_bed_dynamics(
            grid, gsd=gsd, sediment_transport__sediment_supply_imposed_link=testArray
        )


def test_error_gsd_fixed_node():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="bed_surface__grain_size_distribution_fixed_node"
        ".*does not have the same dimensions of the grid's nodes",
    ):
        river_bed_dynamics(
            grid, gsd=gsd, bed_surface__grain_size_distribution_fixed_node=testArray
        )


def test_error_elevation_fixed_node():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="bed_surface__elevation_fixed_node"
        ".*does not have the same dimensions of the grid's nodes",
    ):
        river_bed_dynamics(grid, gsd=gsd, bed_surface__elevation_fixed_node=testArray)


def test_error_bedload_gsd_imposed():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="sediment_transport__bedload_grain_size_distribution_imposed_link"
        ".*does not have the same dimensions of the grid's link",
    ):
        river_bed_dynamics(
            grid,
            gsd=gsd,
            sediment_transport__bedload_grain_size_distribution_imposed_link=testArray,
        )


def test_error_previous_Velocity():
    # Set the topographic elevation array
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.add_zeros("topographic__elevation", at="node")
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    testArray = np.array([0, 1, 2, 3])

    # Check for ValueError when instantiating river_bed_dynamics
    with pytest.raises(
        ValueError,
        match="surface_water__velocity_previous_time_link"
        ".*does not have the same dimensions of the grid's link",
    ):
        river_bed_dynamics(
            grid, gsd=gsd, surface_water__velocity_previous_time_link=testArray
        )


def test_outlet_links_horizontal_right_edge():
    # Set the topographic elevation array
    values = np.full((34, 1), 45.00)
    middle_values = np.arange(24.00, -1.00, -0.75).reshape(-1, 1)
    dem = np.hstack((values, middle_values, middle_values, values)).T
    topographic__elevation = dem.flatten()  # In landlab format
    grid = RasterModelGrid((4, 34), xy_spacing=50)
    grid.at_node["topographic__elevation"] = topographic__elevation
    grid.at_node["topographic__elevation"][101] = -0.8

    gsd = np.array([[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [2, 10], [1, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    grid.set_watershed_boundary_condition(topographic__elevation)

    rbd = river_bed_dynamics(grid, gsd=gsd)
    rbd.run_one_step()

    outlet_links_horizontal_sol = np.array([[165], [166]])
    np.testing.assert_almost_equal(
        rbd._outlet_links_horizontal, outlet_links_horizontal_sol, decimal=1
    )


def test_outlet_links_horizontal_bottom_edge():
    # Set the topographic elevation array
    values = np.full((34, 1), 45.00)
    middle_values = np.arange(24.00, -1.00, -0.75).reshape(-1, 1)
    dem = np.hstack((values, middle_values, middle_values, values))
    topographic__elevation = np.flip(dem, 0).flatten()  # In landlab format
    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.at_node["topographic__elevation"] = topographic__elevation
    grid.at_node["topographic__elevation"][1] = -0.8
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")

    grid.set_watershed_boundary_condition(topographic__elevation)

    rbd = river_bed_dynamics(grid, gsd=gsd)
    rbd.run_one_step()

    outlet_links_horizontal_sol = np.array([[7, 0], [8, 1]])
    np.testing.assert_almost_equal(
        rbd._outlet_links_horizontal, outlet_links_horizontal_sol, decimal=1
    )


def test_r_b_d_approximate_solution():
    # Set the topographic elevation array
    values = np.full((34, 1), 45.00)
    middle_values = np.arange(24.00, -1.00, -0.75).reshape(-1, 1)
    dem = np.hstack((values, middle_values, middle_values, values))
    topographic__elevation = np.flip(dem, 0).flatten()  # In landlab format
    z = topographic__elevation

    # Set Numerical simulation conditions and time control settings
    max_dt = 5
    simulation_max_time = 0.25 * 86400

    # Flow, bed, and upstream simulation conditions
    n = 0.03874  # Manning's n
    upstream_sediment_supply = -0.0087  # bedload rate at inlet

    # Link Id in which sediment supply and discharge enters
    link_inlet = np.array((221, 222))

    # Node Id in Water depth is specified
    node_inlet = np.array((129, 130))

    # Node ID for fixed Nodes
    fixed_nodes_id = np.array((1, 2, 5, 6))

    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.at_node["topographic__elevation"] = topographic__elevation
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    of = OverlandFlow(
        grid,
        h_init=0.001,
        mannings_n=n,
        rainfall_intensity=0.0,
    )
    of._rainfall_intensity = np.zeros_like(z, dtype=float)
    of._rainfall_intensity[
        node_inlet
    ] = 0.02  # Boundary conditions of discharge and flow depth

    # Creates fields and instantiate the RiverbedDynamics component
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity", at="link")
    grid["link"]["surface_water__depth"] = map_mean_of_link_nodes_to_link(
        grid, "surface_water__depth"
    )

    fixed_nodes = np.zeros_like(z)
    fixed_nodes[fixed_nodes_id] = 1
    # sediment_transport__sediment_supply_imposed
    qb = np.full(grid.number_of_links, 0.0)
    qb[link_inlet] = upstream_sediment_supply

    rbd = river_bed_dynamics(
        grid,
        gsd=gsd,
        variable_critical_shear_stress=True,
        outlet_boundary_condition="fixedValue",
        bed_surface__elevation_fixed_node=fixed_nodes,
        sediment_transport__sediment_supply_imposed_link=qb,
    )

    # Set boundaries as closed boundaries, the outlet is set to an open boundary.
    grid.set_watershed_boundary_condition_outlet_id([1, 2], z, 45.0)

    # Initial Conditions
    t0 = 0
    while t0 < 3600:
        of.overland_flow(dt=max_dt)  # Runs overland flow for one time step
        t0 += of.dt

    # Node ID for calculated node elevation
    calculated_nodes_Id = np.arange(128, 136)
    number_columns = grid.number_of_node_columns
    number_rows_calculated_nodes = int(calculated_nodes_Id.shape[0] / number_columns)
    calculated_nodes_Id = np.reshape(
        calculated_nodes_Id, (number_rows_calculated_nodes, number_columns)
    )

    t = 0.0
    while t < simulation_max_time:
        # defines the velocity at previous time
        rbd._surface_water__velocity_previous_time_link = (
            of._grid["link"]["surface_water__discharge"]
            / of._grid["link"]["surface_water__depth"]
        )

        # Runs overland flow for one time step
        of.overland_flow(dt=max_dt)

        # defines the velocity at current time
        grid["link"]["surface_water__velocity"] = (
            grid["link"]["surface_water__discharge"]
            / grid["link"]["surface_water__depth"]
        )

        # Defines the time step used in river_bed_dynamics
        rbd._grid._dt = of.dt

        # Runs river_bed_dynamics for one time step
        rbd.run_one_step()  # Runs riverBedDynamics for one time step

        # Gradient preserving at upstream ghost cells
        dsNodesId = np.array(
            calculated_nodes_Id[0, 1] - np.arange(1, 3) * number_columns
        )
        z = grid["node"]["topographic__elevation"]  # Updated topographic elevation
        bedSlope = (z[dsNodesId[0]] - z[dsNodesId[1]]) / grid.dx

        for i in np.arange(0, calculated_nodes_Id.shape[0]):
            grid["node"]["topographic__elevation"][
                calculated_nodes_Id[i, 1 : number_columns - 1]
            ] = (
                z[calculated_nodes_Id[i, 1 : number_columns - 1] - 2 * number_columns]
                + 2 * grid.dx * bedSlope
            )

        t = t + of.dt

    z = np.reshape(z, dem.shape)[:, 1]
    # x = np.arange(0, 1700, 50, dtype=np.int64)
    z_solution_approx = np.array(
        [
            -0.75000,
            0.00000,
            0.74019,
            1.49101,
            2.24352,
            2.99569,
            3.74731,
            4.49841,
            5.24912,
            5.99954,
            6.74978,
            7.49990,
            8.24997,
            9.00000,
            9.75002,
            10.50000,
            11.25010,
            12.00010,
            12.75030,
            13.50060,
            14.25140,
            15.00330,
            15.75730,
            16.51570,
            17.28220,
            18.06310,
            18.86790,
            19.70970,
            20.60410,
            21.56840,
            22.61440,
            23.78860,
            24.96270,
            26.13690,
        ]
    )

    np.testing.assert_almost_equal(z, z_solution_approx, decimal=1)
