#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 08:42:24 2021

@author: gtucker
"""

from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_raises

from landlab import HexModelGrid, RadialModelGrid, RasterModelGrid
from landlab.components import Flexure, ListricKinematicExtender


def test_hangingwall_nodes():
    """Test the correct identification of hangingwall nodes."""
    grid = RasterModelGrid((3, 7), xy_spacing=2500.0)
    grid.add_zeros("topographic__elevation", at="node")
    extender = ListricKinematicExtender(grid, fault_location=2500.0)

    assert_array_equal(
        extender._hangwall, [2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20]
    )


def test_subsidence_and_horiz_shift():
    """Test that elev subsides then shifts after 2 time steps."""
    grid = RasterModelGrid((3, 7), xy_spacing=2500.0)
    topo = grid.add_zeros("topographic__elevation", at="node")
    extender = ListricKinematicExtender(
        grid, extension_rate=0.01, fault_location=2500.0
    )

    # Run long enough to extend by half a grid cell
    extender.run_one_step(dt=125000.0)
    assert_array_almost_equal(
        topo[7:14],
        [0.0, 0.0, -1404.156819, -910.66907, -590.616478, -383.045648, -248.425118],
    )

    # Now extend another half cell, so cumulative extension is one cell and
    # elevations should get shifted by one cell
    extender.run_one_step(dt=125000.0)
    assert_array_almost_equal(
        topo[7:14],
        [0.0, 0.0, -3514.477461, -2808.313638, -1821.338140, -1181.232956, -766.091296],
    )

    # Another step, and this time the hangingwall edge has moved by one cell,
    # so the first 3 cells in this row should not further subside.
    extender.run_one_step(dt=125000.0)
    assert_array_almost_equal(
        topo[7:14],
        [
            0.0,
            0.0,
            -3514.477461,
            -3718.982708,
            -2411.954617,
            -1564.278603,
            -1014.516414,
        ],
    )


def test_with_hex_grid():
    grid = HexModelGrid((5, 5), node_layout="rect")
    grid.add_zeros("topographic__elevation", at="node")
    ListricKinematicExtender(grid)
    ListricKinematicExtender(grid, fault_location=2.0)

    grid = HexModelGrid((5, 5), node_layout="rect", orientation="vertical")
    grid.add_zeros("topographic__elevation", at="node")
    assert_raises(NotImplementedError, ListricKinematicExtender, grid)


def test_with_flexure():
    """Test integrating with flexure."""
    crust_density = 2700.0  # density of crustal column, kg/m3
    dx = 2500.0  # grid spacing, m
    dt = 125000.0  # time step, y
    upper_crust_base_depth = 10000.0  # m

    grid = RasterModelGrid((3, 7), xy_spacing=dx)
    topo = grid.add_zeros("topographic__elevation", at="node")
    load = grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    thickness = grid.add_zeros("upper_crust_thickness", at="node")
    upper_crust_base = grid.add_zeros("upper_crust_base__elevation", at="node")

    extender = ListricKinematicExtender(
        grid,
        extension_rate=0.01,
        fault_location=2500.0,
        track_crustal_thickness=True,
    )
    flexer = Flexure(grid, eet=5000.0, method="flexure")
    deflection = grid.at_node["lithosphere_surface__elevation_increment"]

    topo[
        grid.x_of_node <= 7500.0
    ] = 1000.0  # this will force thickness to be 1 km greater at left
    upper_crust_base[:] = -upper_crust_base_depth
    thickness[:] = topo - upper_crust_base
    unit_wt = crust_density * flexer.gravity
    load[:] = unit_wt * thickness  # loading pressure

    # Get the initial deflection, which we'll need to calculate total current
    # deflection
    flexer.update()
    init_deflection = deflection.copy()

    # Run extension for half a grid cell. Elevations change, but thickness
    # doesn't, so deflection should not change. We should be able to recover
    # elevation from:
    #
    #   topo = thickness + crust base - (deflection + subsidence)
    #
    extender.run_one_step(dt=dt)
    flexer.update()
    net_deflection = deflection - init_deflection
    assert_array_almost_equal(
        net_deflection[7:14],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    test_topo = thickness + upper_crust_base - (net_deflection + extender._cum_subs)
    assert_array_almost_equal(topo, test_topo)

    # Now extend for another half cell, which should force a shift. The
    # cumulative subsidence will be subtracted from the thickness field,
    # representing thinning as the hangingwall slides to the "right". This
    # will cause net upward isostatic deflection.
    extender.run_one_step(dt=dt)
    load[:] = unit_wt * thickness
    flexer.update()
    net_deflection = deflection - init_deflection
    assert_array_almost_equal(
        thickness[7:14],
        [
            11000.0,
            11000.0,
            8191.686362,  # greatest subsidence: lost nearly 3 km
            9178.66186,
            9818.767044,  # thicker because shifted (only lost <200 m)
            9233.908704,
            9503.149763,
        ],
    )
    assert_array_almost_equal(
        net_deflection[7:14],
        [
            -59.497362,
            -65.176276,
            -69.222531,
            -70.334462,
            -68.608952,
            -64.912352,
            -59.743080,
        ],
    )


def test_error_handling():

    radial_grid = RadialModelGrid(
        n_rings=1, nodes_in_first_ring=8
    )  # , xy_of_center=(0., 0.))
    assert_raises(TypeError, ListricKinematicExtender, radial_grid)

    hex_grid = HexModelGrid((3, 3))
    assert_raises(TypeError, ListricKinematicExtender, hex_grid)

    grid = RasterModelGrid((3, 7))
    grid.add_zeros("topographic__elevation", at="node")
    assert_raises(
        KeyError, ListricKinematicExtender, grid, track_crustal_thickness=True
    )
