#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 10:22:59 2020

@author: gtucker
"""

import numpy as np
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_equal,
    assert_raises,
)

from landlab import HexModelGrid, RadialModelGrid, RasterModelGrid
from landlab.components import TidalFlowCalculator


def test_constant_depth_deeper_than_tidal_amplitude():
    r"""Test velocity calculation under the following conditions:

    r = 1  # tidal range in m
    T = 4.0e4  # tidal period in s
    n = 0.01  # roughness, s/m^1/3
    chi = 1  # scale velocity, m/s
    h = 8    # tidal mean depth, m

    Under these conditions, the key factors are:

    velocity coefficient, Cv = h^(4/3) / n^2 chi = 1.6e5 m/s
    diffusion coefficient, D = h^(7/3) / n^2 chi = 1.28e6 m2/s
    inundation rate (I) = 2 r / T = 5e-5 m/s

    Domain, L: 300 m, 1D

    Analytical solution for water surface height, ebb tide:

        $\eta = (IL^2/D) ((x/L) - (1/2) (x/L)^2)$

    Solution for x-directed velocity. Note that in this case $I$ and $h$ are
    uniform, so

        $u h = I x$  # output = input

        $u(x) = I x / h$

    So with the above parameters, $I$ is 5e-5 m/s. In a grid with 5 columns,
    $x$ is effectively 100, 200, and 300 m, with the adjacent open boundary node
    at 350 m. So the velocity in m/s should be:
        0.000625, 0.00125, 0.001875
    """
    grid = RasterModelGrid((3, 5), xy_spacing=100.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[:] = -8.0
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)

    tfc = TidalFlowCalculator(grid, tidal_period=4.0e4)
    tfc.run_one_step()

    assert_array_almost_equal(
        grid.at_link["ebb_tide_flow__velocity"][10:13], [0.000625, 0.00125, 0.001875]
    )
    assert_array_almost_equal(
        grid.at_link["flood_tide_flow__velocity"][10:13],
        [-0.000625, -0.00125, -0.001875],
    )


def test_constant_depth_deeper_than_tidal_amplitude_alt_grid():
    """Test velocity calculation with different grid orientation."""
    grid = RasterModelGrid((5, 3), xy_spacing=100.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[:] = -8.0
    grid.set_closed_boundaries_at_grid_edges(True, False, True, True)

    tfc = TidalFlowCalculator(grid, tidal_period=4.0e4)
    tfc.run_one_step()

    links_to_test = [8, 13, 18]
    ebb_vel = grid.at_link["ebb_tide_flow__velocity"][links_to_test]
    assert_array_almost_equal(ebb_vel, [0.000625, 0.00125, 0.001875])
    flood_vel = grid.at_link["flood_tide_flow__velocity"][links_to_test]
    assert_array_almost_equal(flood_vel, [-0.000625, -0.00125, -0.001875])


def test_constant_depth_shallower_than_tidal_amplitude():
    r"""Test velocity calculation under the following conditions:

    r = 1  # tidal range in m
    T = 4.0e4  # tidal period in s
    n = 0.01  # roughness, s/m^1/3
    chi = 1  # scale velocity, m/s
    h_m = 0.25    # water depth at mean sea level, m

    Under these conditions, the key factors are:

    inundation rate (I) =

    $I = [r/2 − max(−r/2, min(z, r/2))]/(T/2)$

    $= (0.5 - max(-0.25, min(-0.25, 0.5))) / 20,000$

    $= (0.5 - (-0.25)) / 20,000 = 3.75\times 10^{-5}$ m/s

    Domain, L: 300 m, 1D

    Analytical solution for water surface height, ebb tide:

        $\eta = (IL^2/D) ((x/L) - (1/2) (x/L)^2)$

    Solution for x-directed velocity. Note that in this case $I$ and $h$ are
    uniform, so

        $u h = I x$  # output = input

        $u(x) = I x / h$

    Water depth $h$ is calculated as the average of high-tide and low-tide
    depth. Because low-tide depth is zero,

        $h = (0.25 + r/2) / 2 = 0.375$ m

    In a grid with 5 columns, $x$ is effectively 100, 200, and 300 m, with the
    adjacent open boundary node at 350 m. So the velocity in m/s should be:
        0.01, 0.02, 0.03
    """
    grid = RasterModelGrid((3, 5), xy_spacing=100.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[:] = -0.25
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)

    tfc = TidalFlowCalculator(grid, tidal_period=4.0e4)
    tfc.run_one_step()

    assert_array_almost_equal(
        grid.at_link["ebb_tide_flow__velocity"][10:13], [0.01, 0.02, 0.03]
    )
    assert_array_almost_equal(
        grid.at_link["flood_tide_flow__velocity"][10:13], [-0.01, -0.02, -0.03]
    )


def test_with_hex_grid():
    """Test mass balance with a hex grid.

    The test here is based on a simple mass balance: the computed velocity at
    the open boundaries, when multiplied by the depth and the width of all open
    cell faces, should equal the total inflow or outflow rate, which is the
    inundation rate times the total area.

    The test configuration is a hex grid with 5 rows and a maximum width of 5
    columns. The bottom 3 nodes are open (fixed value) boundaries; the rest are
    closed. The tidal range is 1 meter, the mean depth is 5 meters, and the
    tidal period is 40,000 seconds. Node spacing will be 2 meters.

    Inundation rate = I = tidal range / tidal half period
    = 1 / 20,000 = 5 x 10^-5 m/s

    Area of one cell = (3^0.5 / 2) dx^2 ~ 3.4641

    Width of one face = dx / 3^0.5

    Inundation volume rate = I x cell area x 7 cells
    = ~0.0012124 m3/s

    Outflow volume = velocity at one the edge of any one of the lower active
    links x (solved) depth at that link x width of face x 4 faces. Call the
    velocity-depth product q. Then the predicted q should be:

        q = inundation volume rate / (face width x 4 faces)
          = (7 I (3^0.5 / 2) dx^2) / (4 dx / 3^0.5)
          = (21/8) I dx = (21/4) r dx / T = 0.0002625
    """
    grid = HexModelGrid((5, 3), spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[:] = -5.0

    # Close all boundary nodes except the bottom row
    grid.status_at_node[
        grid.status_at_node != grid.BC_NODE_IS_CORE
    ] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0:3] = grid.BC_NODE_IS_FIXED_VALUE

    tfc = TidalFlowCalculator(grid, tidal_period=4.0e4)
    tfc.run_one_step()

    q = grid.at_link["flood_tide_flow__velocity"] * tfc._water_depth_at_links
    assert_array_almost_equal(q[3:7], [0.0002625, 0.0002625, 0.0002625, 0.0002625])


def test_wrong_grid_type_error():
    grid = RadialModelGrid(2)
    assert_raises(TypeError, TidalFlowCalculator, grid)


def test_getters_and_setters():
    grid = RasterModelGrid((3, 5))
    grid.add_zeros("topographic__elevation", at="node")
    tfc = TidalFlowCalculator(grid)

    tfc.roughness = 0.01
    assert_array_equal(tfc.roughness, 0.01 + np.zeros(grid.number_of_links))

    tfc.tidal_range = 4.0
    assert_equal(tfc.tidal_range, 4.0)
    assert_equal(tfc._tidal_half_range, 2.0)

    tfc.tidal_period = 40000.0
    assert_equal(tfc.tidal_period, 40000.0)
    assert_equal(tfc._tidal_half_period, 20000.0)

    tfc.mean_sea_level = 1.0
    assert_equal(tfc.mean_sea_level, 1.0)
