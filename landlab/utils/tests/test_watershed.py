#!/usr/bin/env python

from landlab import RasterModelGrid
from landlab.components import FlowRouter
from landlab.utils import get_watershed_nodes, get_watershed_outlet
import numpy as np


def test_get_watershed_nodes():

    grid = RasterModelGrid((7, 7), 1)

    z = np.array([
        -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
        -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
        -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
        -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
        -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
        -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
        -9999., -9999., -9999., -9999., -9999., -9999., -9999.])

    grid.at_node['topographic__elevation'] = z

    outlet_id = 2

    grid.set_watershed_boundary_condition_outlet_id(outlet_id, z,
                                                    nodata_value=-9999.)

    fr = FlowRouter(grid)
    fr.run_one_step()

    ws_nodes = get_watershed_nodes(grid, outlet_id)

    # Given the watershed boundary conditions, the number of watershed nodes
    # should be equal to the number of core nodes plus 1 for the outlet node.
    np.testing.assert_equal(len(ws_nodes), grid.number_of_core_nodes + 1)

def test_get_watershed_outlet():

    grid = RasterModelGrid((7, 7), 1)

    z = np.array([
        -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
        -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
        -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
        -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
        -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
        -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
        -9999., -9999., -9999., -9999., -9999., -9999., -9999.])
 
    grid.at_node['topographic__elevation'] = z

    imposed_outlet = 2

    grid.set_watershed_boundary_condition_outlet_id(imposed_outlet, z,
                                                    nodata_value=-9999.)

    fr = FlowRouter(grid)
    fr.run_one_step()

    test_node = 32

    determined_outlet = get_watershed_outlet(grid, test_node)
    np.testing.assert_equal(determined_outlet, imposed_outlet)

    # Create a pit.
    pit_node = 38
    grid.at_node['topographic__elevation'][pit_node] -= 32
    fr.run_one_step()

    pit_outlet = get_watershed_outlet(grid, test_node)
    np.testing.assert_equal(pit_outlet, pit_node)
