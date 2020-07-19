#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 10:22:59 2020

@author: gtucker
"""

from tidal_flow_calculator import TidalFlowCalculator
from numpy.testing import assert_array_almost_equal
from landlab import RasterModelGrid


def test_constant_depth_deeper_than_tidal_amplitude():
    """Test velocity calculation under the following conditions:

    r = 1  # tidal range in m
    T = 4.0e4  # tidal period in s
    n = 0.01  # roughness, s/m^1/3
    chi = 1  # scale velocity, m/s
    h = 8    # tidal mean depth, m

    Under these conditions, the key factors are:

    velocity coefficient, Cv = h^(4/3) / n^2 chi = 1.6e5 m/s
    diffusion coefficient, D = h^(7/3) / n^2 chi = 1.28e6 m2/s
    inundation rate (I) = r / T = 2.5e-5 m/s

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
    z = grid.add_zeros('topographic__elevation', at='node')
    z[:] = -8.0
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)

    tfc = TidalFlowCalculator(grid, tidal_period=4.0e4)
    tfc.run_one_step()

    assert_array_almost_equal(grid.at_link["ebb_tide_flow__velocity"][10:13],
                              [0.000625, 0.00125, 0.001875])
    assert_array_almost_equal(grid.at_link["flood_tide_flow__velocity"][10:13],
                              [-0.000625, -0.00125, -0.001875])


test_constant_depth_deeper_than_tidal_amplitude()