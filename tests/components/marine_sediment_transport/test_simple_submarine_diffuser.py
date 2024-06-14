#!/usr/bin/env python3
"""
Created on Fri Feb 12 17:12:24 2021

@author: gtucker
"""

from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import SimpleSubmarineDiffuser


def test_diffusivity_vs_depth():
    """
    Make sure the diffusivity-vs-depth function is correct.

    This test checks that the diffusivity should equal the shallow water
    diffusivity at the shoreline (when tidal range is zero), and at the wave
    base depth (here 50 m). It should be about half the nominal value at a
    depth of 84.657359 m (1 - ln(1/2) times the wave base depth). It should be
    near-zero for negative depths.

    When the tidal range is nonzero, the value at a depth of zero (shoreline)
    should be half the shallow-water diffusivity.
    """
    grid = RasterModelGrid((3, 3), xy_spacing=40.0)
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[3] = -84.657359
    topo[4] = -50.0
    topo[5] = 50.0

    ssd = SimpleSubmarineDiffuser(
        grid, tidal_range=0.0, wave_base=50.0, shallow_water_diffusivity=100.0
    )
    ssd.calc_diffusion_coef()
    assert_array_almost_equal(grid.at_node["kd"][2:6], [100.0, 50.0, 100.0, 1.0e-20])

    ssd = SimpleSubmarineDiffuser(
        grid, tidal_range=2.0, wave_base=50.0, shallow_water_diffusivity=100.0
    )
    ssd.calc_diffusion_coef()
    assert_array_almost_equal(grid.at_node["kd"][2:6], [50.0, 50.0, 100.0, 1.0e-20])


def test_one_step_shallow():
    r"""
    Test one time step of diffusion in shallow water.

    The test uses a submarine topography above wave base, with a cross-section
    that looks like this:

        o_./'\._o

    where o is an open-boundary node, . and ' indicate a core node, and the
    slope of the / and \ links is 0.1. The other boundary nodes are closed.
    With a diffusivity coefficient of 100 m2/y, the flux along the sloping
    links should be 10 m2/y. Let the node spacing be 40 m. The height of most
    of the nodes will be -10 m, and the central node -6 m. With a time step of
    2 years, the central node should lose 20 m2 to each side, for a total of
    40 m2. This translates into an elevation reduction of 1 m. The nodes on
    either side each gain 40 m2, for an elevation gain of 0.5 m. Therefore the
    3 core nodes should change from elevatios of [-10, -6, -10] to
    [-9.5, -7, -9.5].
    """
    grid = RasterModelGrid((3, 5), xy_spacing=40.0)
    grid.set_closed_boundaries_at_grid_edges(
        right_is_closed=False,
        top_is_closed=True,
        left_is_closed=False,
        bottom_is_closed=True,
    )
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[:] = -10.0
    topo[7] = -6.0

    ssd = SimpleSubmarineDiffuser(
        grid, tidal_range=0.0, wave_base=100.0, shallow_water_diffusivity=100.0
    )
    ssd.run_one_step(dt=2.0)
    assert_array_equal(topo[6:9], [-9.5, -7, -9.5])
    assert_array_equal(
        grid.at_node["sediment_deposit__thickness"][6:9], [0.5, -1.0, 0.5]
    )


def test_one_step_deep():
    """
    Test one time step of diffusion in deep water (below wave base).

    This test is similar to the shallow-water test but operates below the wave
    base. Let the wave base be 50 m. Using the exponential formulation for
    depth-decay of transport, the transport coefficient at depth D is
    K0 exp( -(D - Dwb) / Dwb ), where Dwb is the wave-base depth, and K0 is the
    shallow-water transport coefficient (diffusivity). If D = 85 m, then the
    transport coefficient is about 0.5 K0.
    """
    grid = RasterModelGrid((3, 5), xy_spacing=40.0)
    grid.set_closed_boundaries_at_grid_edges(
        right_is_closed=False,
        top_is_closed=True,
        left_is_closed=False,
        bottom_is_closed=True,
    )
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[:] = -87.0
    topo[7] = -83.0

    ssd = SimpleSubmarineDiffuser(
        grid, tidal_range=0.0, wave_base=50.0, shallow_water_diffusivity=100.0
    )
    ssd.run_one_step(dt=2.0)
    assert_array_almost_equal(topo[6:9], [-86.741574, -83.516851, -86.741574])


def test_depth_function():
    """
    Test the depth-weighting function.

    If tidal range is zero, the weight should be 1 where water depth is >=0,
    and 0 otherwise. If tidal range is >0, the weighting function should be 0.5
    where depth = 0, about 0.8808 where depth equals tidal range, and about
    0.1192 where depth equals minus one tidal range.
    """
    grid = RasterModelGrid((3, 4))
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[5] = -2.0
    topo[6] = 2.0

    ssd = SimpleSubmarineDiffuser(grid, tidal_range=0.0)
    depth = -topo
    df = ssd.depth_function(depth)
    assert_array_equal(df, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    ssd = SimpleSubmarineDiffuser(grid, tidal_range=2.0)
    df = ssd.depth_function(depth)
    assert_array_almost_equal(
        df, [0.5, 0.5, 0.5, 0.5, 0.5, 0.880797, 0.119203, 0.5, 0.5, 0.5, 0.5, 0.5]
    )
