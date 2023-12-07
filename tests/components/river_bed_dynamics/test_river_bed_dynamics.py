"""
Unit tests for landlab.components.RiverBedDynamics.RiverBedDynamics

last updated: 10/12/2023
"""
import numpy as np

from landlab import RasterModelGrid
from landlab.components import OverlandFlow, RiverBedDynamics
from landlab.grid.mappers import map_mean_of_link_nodes_to_link

(_SHAPE, _SPACING, _ORIGIN) = ((5, 5), (100, 100), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(rbd):
    assert rbd.name == "RiverBedDynamics"


def test_input_var_names(rbd):
    assert rbd.input_var_names == (
        "surface_water__depth",
        "surface_water__velocity",
        "topographic__elevation",
    )


def test_output_var_names(rbd):
    assert rbd.output_var_names == ("topographic__elevation",)


def test_optional_var_names(rbd):
    assert rbd.optional_var_names == ()


def test_var_units(rbd):
    assert set(rbd.input_var_names) | set(rbd.output_var_names) | set(
        rbd.optional_var_names
    ) == set(dict(rbd.units).keys())

    assert rbd.var_units("surface_water__depth") == "m"
    assert rbd.var_units("surface_water__velocity") == "m/s"
    assert rbd.var_units("topographic__elevation") == "m"


def test_grid_shape(rbd):
    assert rbd.grid.number_of_node_rows == _SHAPE[0]
    assert rbd.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(rbd):
    assert rbd.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(rbd):
    assert rbd.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_median_size(rbd):
    assert rbd._bed_surf__median_size_node[20] == 16


def test_geometric_mean_size_node(rbd):
    np.testing.assert_almost_equal(rbd._bed_surf__geo_std_size_node, 2.606, decimal=2)


def test_sand_fraction_node(rbd):
    np.testing.assert_almost_equal(rbd._bed_surf__sand_fract_node[20], 0.01, decimal=3)


def test_velocity_previous_time(rbd):
    np.testing.assert_almost_equal(
        rbd._surface_water__velocity_prev_time_link[20], 10, decimal=3
    )


def test_rbd_approximate_solution():
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

    # Set boundaries as closed boundaries, the outlet is set to an open boundary.
    grid.set_watershed_boundary_condition_outlet_id([1, 2], z, 45.0)

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

    rbd = RiverBedDynamics(
        grid,
        gsd=gsd,
        variable_critical_shear_stress=True,
        outlet_boundary_condition="fixedValue",
        bed_surf__elev_fix_node=fixed_nodes,
        sed_transp__bedload_rate_fix_link=qb,
    )

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
        rbd._surface_water__velocity_prev_time_link = (
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

        # Defines the time step used in RiverBedDynamics
        rbd._grid._dt = of.dt

        # Runs RiverBedDynamics for one time step
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
            -0.7500000,
            0.00000000,
            0.74017691,
            1.49102685,
            2.24353750,
            2.99571190,
            3.74732126,
            4.49841756,
            5.24911398,
            5.99952878,
            6.74976148,
            7.49988502,
            8.24994754,
            8.99997843,
            9.74999517,
            10.50000965,
            11.25003498,
            12.00009545,
            12.75024571,
            13.50061034,
            14.25146400,
            15.00338512,
            15.75753064,
            16.51608576,
            17.28291946,
            18.06439036,
            18.87007451,
            19.71299066,
            20.60881050,
            21.57444743,
            22.62152078,
            23.79589244,
            24.97026410,
            26.14463576,
        ]
    )

    np.testing.assert_almost_equal(z, z_solution_approx, decimal=1)
