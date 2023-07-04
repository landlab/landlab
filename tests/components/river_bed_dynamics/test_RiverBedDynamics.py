"""
Unit tests for landlab.components.river_bed_dynamics.RiverBedDynamics

last updated: 06/09/2023
"""
import numpy as np

from landlab import RasterModelGrid
from landlab.components import OverlandFlow, RiverBedDynamics
from landlab.grid.mappers import map_mean_of_link_nodes_to_link

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(r_b_d):
    assert r_b_d.name == "RiverBedDynamics"


def test_input_var_names(r_b_d):
    assert r_b_d.input_var_names == (
        "bed_surface__grain_size_distribution_location",
        "surface_water__depth",
        "surface_water__velocity",
        "surface_water__velocity_previous_time",
        "topographic__elevation",
    )


def test_output_var_names(r_b_d):
    assert r_b_d.output_var_names == (
        "bed_surface__geometric_mean_size",
        "bed_surface__geometric_standard_deviation_size",
        "bed_surface__grain_size_distribution",
        "bed_surface__median_size",
        "bed_surface__sand_fraction",
        "sediment_transport__bedload_grain_size_distribution",
        "sediment_transport__bedload_rate",
        "sediment_transport__net_bedload",
        "surface_water__shear_stress",
    )


def test_optional_var_names(r_b_d):
    assert r_b_d.optional_var_names == (
        "bed_surface__elevation_fixed",
        "bed_surface__grain_size_distribution_fixed",
        "sediment_transport__bedload_grain_size_distribution_imposed",
        "sediment_transport__sediment_supply_imposed",
    )


def test_var_units(r_b_d):
    assert set(r_b_d.input_var_names) | set(r_b_d.output_var_names) | set(
        r_b_d.optional_var_names
    ) == set(dict(r_b_d.units).keys())

    assert r_b_d.var_units("bed_surface__elevation_fixed") == "-"
    assert r_b_d.var_units("bed_surface__geometric_mean_size") == "mm"
    assert r_b_d.var_units("bed_surface__geometric_standard_deviation_size") == "mm"
    assert r_b_d.var_units("bed_surface__grain_size_distribution") == "-"
    assert r_b_d.var_units("bed_surface__grain_size_distribution_fixed") == "-"
    assert r_b_d.var_units("bed_surface__grain_size_distribution_location") == "-"
    assert r_b_d.var_units("bed_surface__median_size") == "mm"
    assert r_b_d.var_units("bed_surface__sand_fraction") == "-"
    assert r_b_d.var_units("sediment_transport__bedload_grain_size_distribution") == "-"
    assert (
        r_b_d.var_units("sediment_transport__bedload_grain_size_distribution_imposed")
        == "-"
    )
    assert r_b_d.var_units("sediment_transport__bedload_rate") == "m^2/s"
    assert r_b_d.var_units("sediment_transport__net_bedload") == "m^3/s"
    assert r_b_d.var_units("sediment_transport__sediment_supply_imposed") == "m^2/s"
    assert r_b_d.var_units("surface_water__depth") == "m"
    assert r_b_d.var_units("surface_water__velocity") == "m/s"
    assert r_b_d.var_units("surface_water__velocity_previous_time") == "m/s"
    assert r_b_d.var_units("surface_water__shear_stress") == "Pa"
    assert r_b_d.var_units("topographic__elevation") == "m"


def test_grid_shape(r_b_d):
    assert r_b_d.grid.number_of_node_rows == _SHAPE[0]
    assert r_b_d.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(r_b_d):
    assert r_b_d.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(r_b_d):
    assert r_b_d.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_r_b_d_approximate_solution():
    # Set the topographic elevation array
    values = np.full((34, 1), 45.00)
    middle_values = np.arange(24.00, -1.00, -0.75).reshape(-1, 1)
    dem = np.hstack((values, middle_values, middle_values, values))
    topographic__elevation = np.flip(dem, 0).flatten()  # In landlab format

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
    fixed_nodes_Id = np.array((1, 2, 5, 6))

    grid = RasterModelGrid((34, 4), xy_spacing=50)
    grid.at_node["topographic__elevation"] = topographic__elevation
    gsd = np.array([[51, 100], [50, 50], [49, 0]])

    # Creates fields and instantiate the OverlandFlow component
    grid.add_zeros("surface_water__depth", at="node")
    of = OverlandFlow(
        grid, dt_max=max_dt, h_init=0.001, mannings_n=n, use_user_defined_time_step=True
    )

    # Creates fields and instantiate the RiverbedDynamics component
    z = topographic__elevation
    grid["node"]["bed_surface__grain_size_distribution_location"] = np.zeros_like(z)
    grid.add_zeros("surface_water__velocity", at="node")
    grid.add_zeros("surface_water__velocity_previous_time", at="node")
    grid["link"]["surface_water__depth"] = map_mean_of_link_nodes_to_link(
        grid, "surface_water__depth"
    )
    grid["link"]["surface_water__velocity"] = map_mean_of_link_nodes_to_link(
        grid, "surface_water__velocity"
    )
    grid["link"][
        "surface_water__velocity_previous_time"
    ] = map_mean_of_link_nodes_to_link(grid, "surface_water__velocity_previous_time")
    rbd = RiverBedDynamics(
        grid,
        gsd=gsd,
        variable_critical_shear_stress=True,
        outlet_boundary_condition="fixedValue",
    )

    # Set boundaries as closed boundaries, the outlet is set to an open boundary.
    grid.set_watershed_boundary_condition_outlet_id([1, 2], z, 45.0)

    # Creates the fixed nodes information
    fixed_nodes = np.zeros_like(
        z
    )  # fixed_nodes defines as 1 if a node is fixed or 0 if it can vary in elevation
    fixed_nodes[fixed_nodes_Id] = 1
    grid["node"][
        "bed_surface__elevation_fixed"
    ] = fixed_nodes  # Assigns fixed locations to landlab grid

    # Create bed and flow initial condition
    # Create bed and flow initial condition
    grid["link"]["surface_water__discharge"][
        link_inlet
    ] = 50  # Flow discharge in m3/s/m
    grid["node"]["surface_water__depth"][node_inlet] = 0.45  # Flow depth in nodes in m
    grid["link"]["surface_water__depth"][link_inlet] = 0.45  # Flow depth in links in m
    grid["link"]["sediment_transport__sediment_supply_imposed"][
        link_inlet
    ] = upstream_sediment_supply

    # Node ID for calculated node elevation
    calculated_nodes_Id = np.arange(128, 136)
    number_columns = grid.number_of_node_columns
    number_rows_calculated_nodes = int(calculated_nodes_Id.shape[0] / number_columns)
    calculated_nodes_Id = np.reshape(
        calculated_nodes_Id, (number_rows_calculated_nodes, number_columns)
    )

    t = 0.0
    while t < simulation_max_time:
        # Boundary conditions of dishcarge and flow depth
        grid["link"]["surface_water__discharge"][
            link_inlet
        ] = -1  # Flow discharge in m3/s/m
        grid["node"]["surface_water__depth"][
            node_inlet
        ] = 0.45  # Flow depth in nodes in m
        grid["link"]["surface_water__depth"][
            link_inlet
        ] = 0.45  # Flow depth in links in m

        # Velocity at previous time
        grid["link"]["surface_water__velocity_previous_time"] = of._q / of._h_links

        of.overland_flow()  # Runs overland flow for one time step

        # Velocity at current time
        grid["link"]["surface_water__velocity"] = of._q / of._h_links

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
    x = np.arange(0, 1700, 50, dtype=np.int64)
    z_solution_approx = (
        2.066e-12 * x**4
        - 5.579e-09 * x**3
        + 4.714e-06 * x**2
        + 1.370e-02 * x
        - 6.830e-01
    )
    np.testing.assert_almost_equal(z, z_solution_approx, decimal=1)
