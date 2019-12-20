#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:26:31 2019

@author: G Tucker, D Litwin
"""

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import FlowAccumulator, GroundwaterDupuitPercolator


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
    for i in range(100):
        gdp.run_one_step(1e3)

    assert_equal(np.round(gdp._thickness[4], 5), 0.001)


def test_simple_surface_leakage():
    """ test a one-node steady simulation for surface leakage.

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

    for i in range(1000):
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
    gdp = GroundwaterDupuitPercolator(
        rg, recharge_rate=1.0e-8, hydraulic_conductivity=0.01, courant_coefficient=0.01
    )
    for i in range(10):
        gdp.run_with_adaptive_time_step_solver(1e4)

    assert_equal(np.round(gdp._thickness[4], 5), 0.001)


def test_conservation_of_mass_adaptive_dt():
    """ test conservation of mass in a sloping aquifer.

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
    for i in range(500):
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
    """ test that water table is symmetric under constant recharge

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
    for i in range(1000):
        gdp.run_one_step(1e3)

    tc = hmg.at_node["aquifer__thickness"]
    assert_almost_equal(tc[5], tc[31])  # SW-NE
    assert_almost_equal(tc[29], tc[7])  # NW-SE
    assert_almost_equal(tc[16], tc[20])  # W-E


def test_wt_above_surface_standard_run_step():
    """ test that water tables above the topogrpahic elevation are
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
    elev = grid.add_ones("node", "topographic__elevation")
    wt = grid.add_zeros("node", "water_table__elevation")
    _ = grid.add_zeros("aquifer_base__elevation", at="node")
    wt[:] = elev + 1

    # initialize the groundwater model
    gdp = GroundwaterDupuitPercolator(grid, recharge_rate=0.0)
    gdp.run_one_step(1)
    assert_equal(wt[4], 1)


def test_wt_above_surface_adaptive_run_step():
    grid = RasterModelGrid((3, 3))
    grid.set_closed_boundaries_at_grid_edges(True, True, True, False)
    elev = grid.add_ones("node", "topographic__elevation")
    wt = grid.add_zeros("node", "water_table__elevation")
    _ = grid.add_zeros("aquifer_base__elevation", at="node")

    wt[:] = elev + 1

    # initialize the groundwater model
    gdp = GroundwaterDupuitPercolator(grid, recharge_rate=0.0)
    gdp.run_with_adaptive_time_step_solver(1)
    assert_equal(wt[4], 1)
