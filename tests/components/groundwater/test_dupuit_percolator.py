#!/usr/bin/env python3
"""
Created on Tue Jun  4 16:26:31 2019

@author: G Tucker, D Litwin
"""

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import GroundwaterDupuitPercolator
from landlab.grid.mappers import map_mean_of_link_nodes_to_link


def test_simple_water_table():
    """Test a one-node steady simulation.

    Notes
    -----
    The analytical solution for one interior cell with one open boundary is as
    follows. The incoming recharge must equal the outgoing discharge. The
    incoming recharge is R, and the surface area is 1 m2, so the incoming
    volume per unit time is R m/s x 1 m x 1 m. The outgoing discharge is equal
    to the conductivity, K (m/s), times the thickness at the boundary, Hb,
    times the hydraulic gradient, which in this case is (H - 0) / dx = H.
    The model uses the upwind thickness, which in this case is H. Therefore:

        K H^2 = R, or

        H = sqrt( R / K ).

    With R = 10^-8 m/s and K = 10^-2 m/s, we should have

        H = 0.001 m.
    """
    boundaries = {"top": "closed", "left": "closed", "bottom": "closed"}
    rg = RasterModelGrid((3, 3), bc=boundaries)
    rg.add_zeros("aquifer_base__elevation", at="node")
    rg.add_ones("topographic__elevation", at="node")
    gdp = GroundwaterDupuitPercolator(
        rg, recharge_rate=1.0e-8, hydraulic_conductivity=0.01
    )
    for _ in range(100):
        gdp.run_one_step(1e3)

    assert_equal(np.round(gdp._thickness[4], 5), 0.001)


def test_simple_surface_leakage():
    """test a one-node steady simulation for surface leakage.

    Notes
    ----
    This test demonstrates that at steady state when no flow is
    allowed to leave the domain through the subsurface, the surface
    water flux is equal to the recharge flux.

    """
    grid = RasterModelGrid((3, 3), xy_spacing=1.0)
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    grid.add_zeros("aquifer_base__elevation", at="node")
    grid.add_ones("topographic__elevation", at="node")
    gdp = GroundwaterDupuitPercolator(grid, recharge_rate=1.0e-6)

    for _ in range(1000):
        gdp.run_one_step(1e3)

    assert_almost_equal(gdp._qs[4], 1e-6)


def test_simple_water_table_adaptive_dt():
    """Test a one-node steady simulation.

    Notes
    -----
    This test demonstrates the same simple water table as
    test_simple_water_table, but with the run_with_adaptive_time_step_solver
    method.


    """
    boundaries = {"top": "closed", "left": "closed", "bottom": "closed"}
    rg = RasterModelGrid((3, 3), bc=boundaries)
    rg.add_zeros("aquifer_base__elevation", at="node")
    rg.add_ones("topographic__elevation", at="node")
    rg.add_zeros("water_table__elevation", at="node")
    rg.at_node["water_table__elevation"][rg.core_nodes] += 1e-10
    gdp = GroundwaterDupuitPercolator(
        rg,
        recharge_rate=1.0e-8,
        hydraulic_conductivity=0.01,
    )
    for _ in range(10):
        gdp.run_with_adaptive_time_step_solver(1e4)

    assert_equal(np.round(gdp._thickness[4], 5), 0.001)


def test_conservation_of_mass_adaptive_dt():
    """test conservation of mass in a sloping aquifer.

    Notes
    ----
    This test demonstrates conservation of mass from a sloping aquifer
    with constant recharge using the flux and storage calculation methods.
    Note that there is a slight loss of mass as the adaptive timestep solver
    does not calculate fluxes in the intermediate substeps, so accuracy in
    this case is only to ~3 significant digits.
    """

    grid = RasterModelGrid((3, 10), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(True, True, False, True)
    elev = grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("aquifer_base__elevation", at="node")

    elev[:] = grid.x_of_node / 100 + 1
    wt = grid.add_zeros("water_table__elevation", at="node")
    wt[:] = elev

    # initialize the groundwater model
    gdp = GroundwaterDupuitPercolator(
        grid,
        hydraulic_conductivity=0.0005,
        recharge_rate=1e-7,
        courant_coefficient=0.01,
    )
    fa = FlowAccumulator(grid, runoff_rate="surface_water__specific_discharge")

    # initialize fluxes we will record
    recharge_flux = 0
    gw_flux = 0
    sw_flux = 0
    storage_0 = gdp.calc_total_storage()

    dt = 1e4
    for _ in range(500):
        gdp.run_with_adaptive_time_step_solver(dt)
        fa.run_one_step()

        recharge_flux += gdp.calc_recharge_flux_in() * dt
        gw_flux += gdp.calc_gw_flux_out() * dt
        sw_flux += gdp.calc_sw_flux_out() * dt
    storage = gdp.calc_total_storage()

    assert_almost_equal(
        (gw_flux + sw_flux + storage - storage_0) / recharge_flux, 1.0, decimal=3
    )


def test_symmetry_of_solution():
    """test that water table is symmetric under constant recharge

    Notes:
    ----
    Under constant recharge with radially symmetric aquifer base elevation,
    the model should produce a water table that is radially symmetric. This test
    demonstrates this is the case where topographic elevation and aquifer
    base elevation are radially symmetric parabolas on a hexagonal grid.

    """
    hmg = HexModelGrid(shape=(7, 4), spacing=10.0)
    x = hmg.x_of_node
    y = hmg.y_of_node
    elev = hmg.add_zeros("topographic__elevation", at="node")
    elev[:] = 1e-3 * (x * (max(x) - x) + y * (max(y) - y)) + 2
    base = hmg.add_zeros("aquifer_base__elevation", at="node")
    base[:] = elev - 2
    wt = hmg.add_zeros("water_table__elevation", at="node")
    wt[:] = elev
    wt[hmg.open_boundary_nodes] = 0.0

    gdp = GroundwaterDupuitPercolator(
        hmg, recharge_rate=1e-7, hydraulic_conductivity=1e-4
    )
    for _ in range(1000):
        gdp.run_one_step(1e3)

    tc = hmg.at_node["aquifer__thickness"]
    assert_almost_equal(tc[5], tc[31])  # SW-NE
    assert_almost_equal(tc[29], tc[7])  # NW-SE
    assert_almost_equal(tc[16], tc[20])  # W-E


def test_wt_above_surface_standard_run_step():
    """test that water tables above the topogrpahic elevation are
    set to the topographic elevation.

    Notes:
    ----
    Water tables above the land surface represent a non-physical condition.
    The GroundwaterDupuitPercolator will not produce this state when it
    is the only component operating on water table elevation or topogrpahic
    elevation, however, when combined with components that do, this may occur.
    If the water table is above the ground surface at a node, it is set to the ground
    surface elevation at that node.

    """

    grid = RasterModelGrid((3, 3))
    grid.set_closed_boundaries_at_grid_edges(True, True, True, False)
    wt = grid.add_ones("node", "water_table__elevation")
    _ = grid.add_ones("node", "topographic__elevation")
    _ = grid.add_zeros("node", "aquifer_base__elevation")

    # initialize the groundwater model
    gdp = GroundwaterDupuitPercolator(grid, recharge_rate=0.0)
    gdp.run_one_step(1)
    assert_equal(wt[4], 1)


def test_wt_above_surface_adaptive_run_step():
    grid = RasterModelGrid((3, 3))
    grid.set_closed_boundaries_at_grid_edges(True, True, True, False)
    wt = grid.add_ones("node", "water_table__elevation")
    _ = grid.add_ones("node", "topographic__elevation")
    _ = grid.add_zeros("node", "aquifer_base__elevation")

    # initialize the groundwater model
    gdp = GroundwaterDupuitPercolator(grid, recharge_rate=0.0)
    gdp.run_with_adaptive_time_step_solver(1)
    assert_equal(wt[4], 1)


def test_inactive_interior_node():
    """
    Test that component returns correct values for recharge flux and
    storage when an interior node is INACTIVE

    Notes:
    ----
    When an interior node is inactive, the number of core nodes is not
    equal to the number of cells. This test confirms that the methods
    to calculate recharge flux and active storage acknowledge this difference.

    """

    mg = RasterModelGrid((4, 4), xy_spacing=1.0)
    mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    elev = mg.add_zeros("node", "topographic__elevation")
    elev[:] = 1
    base = mg.add_zeros("node", "aquifer_base__elevation")
    base[:] = 0
    wt = mg.add_zeros("node", "water_table__elevation")
    wt[:] = 1

    gdp = GroundwaterDupuitPercolator(mg)
    assert_almost_equal(gdp.calc_recharge_flux_in(), 3e-8)
    assert_almost_equal(gdp.calc_total_storage(), 0.6)


def test_k_func():
    """
    Test the use of a function to specify how hydraulic conductivity changes
    with water table position.

    Note:
    ----
    Test that component keeps hydraulic conductivity at default value when
    k_func is None. Then test that a simple function correctly sets the
    hydraulic conductivity value after run_one_step and run_with_adaptive_time_step_solver.

    """

    # initialize model grid
    mg = RasterModelGrid((4, 4), xy_spacing=1.0)
    elev = mg.add_zeros("node", "topographic__elevation")
    elev[:] = 1
    base = mg.add_zeros("node", "aquifer_base__elevation")
    base[:] = 0
    wt = mg.add_zeros("node", "water_table__elevation")
    wt[:] = 0.5

    # initialize model without giving k_func
    gdp = GroundwaterDupuitPercolator(mg)

    # run model and assert that K hasn't changed from the default value
    gdp.run_one_step(0)
    assert np.equal(0.001, gdp.K).all()

    gdp.run_with_adaptive_time_step_solver(0)
    assert np.equal(0.001, gdp.K).all()

    # create a simple k_func, where hydraulic conductivity varies linearly
    # with depth, from Ks at surface to 0 at aquifer base
    def k_func_test(grid, Ks=0.01):
        h = grid.at_node["aquifer__thickness"]
        b = (
            grid.at_node["topographic__elevation"]
            - grid.at_node["aquifer_base__elevation"]
        )
        blink = map_mean_of_link_nodes_to_link(grid, b)
        hlink = map_mean_of_link_nodes_to_link(grid, h)

        return (hlink / blink) * Ks

    # initialize model with given k_func
    gdp1 = GroundwaterDupuitPercolator(mg, hydraulic_conductivity=k_func_test)

    # run model and assert that K has been updated correctly
    gdp1.run_one_step(0)
    assert np.equal(0.005, gdp1.K).all()

    gdp1.run_with_adaptive_time_step_solver(0)
    assert np.equal(0.005, gdp1.K).all()


def test_callback_func():
    """
    Test the use of a callback function to return the storage and
    substep durations while using the run_with_adaptive_time_step_solver
    method.

    Notes:
    ----
    Two tests here: make sure that the substeps sum to the global timestep,
    and make sure that when recharge is 0.0, the total storage does not
    increase during any of the substeps. See component documentation for
    more detail on arguments for the callback_fun.
    """

    # make a function that writes storage and substep duration to
    # externally defined lists
    storage_subdt = []
    subdt = []
    all_n = []

    def test_fun(grid, recharge, dt, n=0.2):
        cores = grid.core_nodes
        h = grid.at_node["aquifer__thickness"]
        area = grid.cell_area_at_node
        storage = np.sum(n * h[cores] * area[cores])

        storage_subdt.append(storage)
        subdt.append(dt)
        all_n.append(n)

    # initialize grid
    grid = RasterModelGrid((3, 3))
    grid.set_closed_boundaries_at_grid_edges(True, True, False, True)
    elev = grid.add_ones("topographic__elevation", at="node")
    elev[3] = 0.1
    grid.add_zeros("aquifer_base__elevation", at="node")
    wt = grid.add_zeros("water_table__elevation", at="node")
    wt[:] = elev

    # initialize groundwater model
    gdp = GroundwaterDupuitPercolator(
        grid,
        recharge_rate=0.0,
        hydraulic_conductivity=0.0001,
        callback_fun=test_fun,
        n=0.1,
    )

    # run groundawter model
    gdp.run_with_adaptive_time_step_solver(1e5)

    # assert that the water table does not increase during substeps
    assert (np.diff(storage_subdt) <= 0.0).all()

    # assert that substeps sum to the global timestep
    assert_almost_equal(1e5, sum(subdt))

    assert all(x == 0.1 for x in all_n)
