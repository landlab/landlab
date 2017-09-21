#!/usr/bin/env python

from landlab import RasterModelGrid
from landlab.components import FlowRouter
from landlab.utils.watershed import get_watershed_ids
import numpy as np

def test_get_watershed_array():
    
    grid = RasterModelGrid((5, 3), 1)
    
    z = np.array([
            -9999., -9999., -9999.,
            -9999.,     1., -9999.,
            -9999.,     5., -9999.,
            -9999.,    32., -9999.,
            -9999., -9999., -9999.])
    
    grid.at_node['topographic__elevation'] = z
    
    outlet_id = 1
    
    grid.set_watershed_boundary_condition_outlet_id(outlet_id, z,
                                                    nodata_value=-9999.)
    
    fr = FlowRouter(grid)
    fr.run_one_step()
    
    ws_ids = get_watershed_ids(grid, outlet_id)
    
    # Given the watershed boundary conditions, the number of watershed ids
    # should be equal to the number of core nodes plus 1 for the outlet node.
    np.testing.assert_equal(len(ws_ids), grid.number_of_core_nodes + 1)
    