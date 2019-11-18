#!/usr/bin/env python
import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.utils import (
    get_watershed_mask,
    get_watershed_masks_with_area_threshold,
    get_watershed_nodes,
    get_watershed_outlet,
)


def test_get_watershed_nodes():
    grid = RasterModelGrid((7, 7))

    z = np.array(
        [
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            26.0,
            0.0,
            30.0,
            32.0,
            34.0,
            -9999.0,
            -9999.0,
            28.0,
            1.0,
            25.0,
            28.0,
            32.0,
            -9999.0,
            -9999.0,
            30.0,
            3.0,
            3.0,
            11.0,
            34.0,
            -9999.0,
            -9999.0,
            32.0,
            11.0,
            25.0,
            18.0,
            38.0,
            -9999.0,
            -9999.0,
            34.0,
            32.0,
            34.0,
            36.0,
            40.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
        ]
    )

    grid.at_node["topographic__elevation"] = z

    outlet_id = 2

    grid.set_watershed_boundary_condition_outlet_id(outlet_id, z, nodata_value=-9999.0)

    fr = FlowAccumulator(grid, flow_director="D8")
    fr.run_one_step()

    ws_nodes = get_watershed_nodes(grid, outlet_id)

    # Given the watershed boundary conditions, the number of watershed nodes
    # should be equal to the number of core nodes plus 1 for the outlet node.
    np.testing.assert_equal(len(ws_nodes), grid.number_of_core_nodes + 1)


def test_get_watershed_masks_with_area_threshold():
    rmg = RasterModelGrid((7, 7), xy_spacing=200)

    z = np.array(
        [
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            26.0,
            0.0,
            26.0,
            30.0,
            34.0,
            -9999.0,
            -9999.0,
            28.0,
            1.0,
            28.0,
            5.0,
            32.0,
            -9999.0,
            -9999.0,
            30.0,
            3.0,
            30.0,
            10.0,
            34.0,
            -9999.0,
            -9999.0,
            32.0,
            11.0,
            32.0,
            15.0,
            38.0,
            -9999.0,
            -9999.0,
            34.0,
            32.0,
            34.0,
            36.0,
            40.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
        ]
    )

    rmg.at_node["topographic__elevation"] = z

    rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    # Route flow.
    fr = FlowAccumulator(rmg, flow_director="D8")
    fr.run_one_step()

    # Get the masks of watersheds greater than or equal to 80,000
    # square-meters.
    critical_area = 80000
    mask = get_watershed_masks_with_area_threshold(rmg, critical_area)

    # Assert that mask null nodes have a drainage area below critical area.
    null_nodes = np.where(mask == -1)[0]
    A = rmg.at_node["drainage_area"][null_nodes]
    below_critical_area_nodes = A < critical_area
    trues = np.ones(len(A), dtype=bool)

    np.testing.assert_equal(below_critical_area_nodes, trues)


def test_get_watershed_outlet():
    grid = RasterModelGrid((7, 7))

    z = np.array(
        [
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            26.0,
            0.0,
            30.0,
            32.0,
            34.0,
            -9999.0,
            -9999.0,
            28.0,
            1.0,
            25.0,
            28.0,
            32.0,
            -9999.0,
            -9999.0,
            30.0,
            3.0,
            3.0,
            11.0,
            34.0,
            -9999.0,
            -9999.0,
            32.0,
            11.0,
            25.0,
            18.0,
            38.0,
            -9999.0,
            -9999.0,
            34.0,
            32.0,
            34.0,
            36.0,
            40.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
            -9999.0,
        ]
    )

    grid.at_node["topographic__elevation"] = z

    imposed_outlet = 2

    grid.set_watershed_boundary_condition_outlet_id(
        imposed_outlet, z, nodata_value=-9999.0
    )

    fr = FlowAccumulator(grid, flow_director="D8")
    fr.run_one_step()

    test_node = 32

    determined_outlet = get_watershed_outlet(grid, test_node)
    np.testing.assert_equal(determined_outlet, imposed_outlet)

    # Create a pit.
    pit_node = 38
    grid.at_node["topographic__elevation"][pit_node] -= 32
    fr.run_one_step()

    pit_outlet = get_watershed_outlet(grid, test_node)
    np.testing.assert_equal(pit_outlet, pit_node)


def test_route_to_multiple_error_raised_watershed_outlet():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        get_watershed_outlet(mg, 10)


def test_route_to_multiple_error_raised_watershed_mask():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        get_watershed_mask(mg, 10)
