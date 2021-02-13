#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 17:12:24 2021

@author: gtucker
"""

from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab import RasterModelGrid
from landlab.components import SimpleSubmarineDiffuser


def test_depth_function():
    """
    Test the depth-weighting function.
    
    If tidal range is zero, the weight should be 1 where water depth is >=0,
    and 0 otherwise. If tidal range is >0, the weighting function should be 0.5
    where depth = 0, about 0.8808 where depth equals tidal range, and about
    0.1192 where depth equals minus one tidal range.
    """
    grid = RasterModelGrid((3, 4))
    topo = grid.add_zeros('topographic__elevation', at='node')
    topo[5] = -2.0
    topo[6] = 2.0

    ssd = SimpleSubmarineDiffuser(grid, tidal_range=0.0)
    depth = -topo
    df = ssd.depth_function(depth)
    assert_array_equal(df, [1.0, 1.0, 1.0, 1.0,
                            1.0, 1.0, 0.0, 1.0,
                            1.0, 1.0, 1.0, 1.0])

    ssd = SimpleSubmarineDiffuser(grid, tidal_range=2.0)
    df = ssd.depth_function(depth)
    assert_array_almost_equal(df, [0.5, 0.5,      0.5,      0.5,
                                   0.5, 0.880797, 0.119203, 0.5,
                                   0.5, 0.5,      0.5,      0.5])


if __name__ == '__main__':
    test_depth_function()

    
    