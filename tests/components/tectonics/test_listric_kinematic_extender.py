#!/usr/bin/env python3
"""
Created on Fri Mar  5 08:42:24 2021

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_raises

from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab.components import ListricKinematicExtender


def test_hangingwall_motion():
    dx = 10.0  # node spacing, m
    fault_x0 = 100.0
    grid = RasterModelGrid((3, 130), xy_spacing=dx)
    elev = grid.add_zeros("topographic__elevation", at="node")
    extender = ListricKinematicExtender(grid, fault_x0=fault_x0, fault_strike=90.0)

    # Verify that node 240 is at xf=1000 and node 190 is at xf=500
    assert_equal(grid.x_of_node[240] - fault_x0, 1000.0)
    assert_equal(grid.x_of_node[190] - fault_x0, 500.0)

    H_init_500 = grid.at_node["hangingwall__thickness"][190]

    dt = 0.2 * dx / 0.001  # time-step duration, y
    nsteps = int(500000.0 / dt)
    print(nsteps)
    for _ in range(nsteps):
        extender.run_one_step(dt)
    assert_almost_equal(
        grid.at_node["hangingwall__thickness"][240], H_init_500, decimal=4
    )
    assert_almost_equal(elev[240], -760.7638, decimal=4)


def test_error_handling():
    radial_grid = RadialModelGrid(n_rings=1, nodes_in_first_ring=8)
    assert_raises(TypeError, ListricKinematicExtender, radial_grid)


def test_preexisting_fields():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("advection__velocity", at="link")
    ListricKinematicExtender(
        grid,
        extension_rate_x=0.0,
        extension_rate_y=1.0,
        fields_to_advect=["hangingwall__thickness"],
    )
    assert_array_almost_equal(
        grid.at_link["advection__velocity"][grid.horizontal_links],
        np.zeros(len(grid.horizontal_links)),
    )
    assert_array_almost_equal(
        grid.at_link["advection__velocity"][grid.vertical_links],
        np.ones(len(grid.vertical_links)),
    )
