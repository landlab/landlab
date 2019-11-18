"""Test the flow router component.

This just tests the functionality of the router in toto - no attempt is made
to test the submodules.
Sinks are tested as part of the lake_mapper.

@author: dejh
"""
# Created on Thurs Nov 12, 2015
import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from six.moves import range

from landlab import CLOSED_BOUNDARY, RadialModelGrid, RasterModelGrid
from landlab.components.flow_routing import FlowRouter

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_deprecation_raised():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    with pytest.deprecated_call():
        FlowRouter(mg)


def test_check_fields(dans_grid1):
    """Check to make sure the right fields have been created."""
    with pytest.deprecated_call():
        FlowRouter(dans_grid1.mg)
    assert_array_equal(dans_grid1.z, dans_grid1.mg.at_node["topographic__elevation"])
    assert_array_equal(np.zeros(25), dans_grid1.mg.at_node["drainage_area"])
    assert_array_equal(np.ones(25), dans_grid1.mg.at_node["water__unit_flux_in"])
    with pytest.deprecated_call():
        FlowRouter(dans_grid1.mg, dans_grid1.infile)
    assert_array_equal(np.full(25, 2.0), dans_grid1.mg.at_node["water__unit_flux_in"])


def test_check_field_input(dans_grid1):
    """Check we can successfully pass water__discharge_in."""
    dans_grid1.mg.add_field(
        "node", "water__unit_flux_in", np.full(25, 3.0), units="m**3/s"
    )
    with pytest.deprecated_call():
        FlowRouter(dans_grid1.mg)
    assert_array_equal(np.full(25, 3.0), dans_grid1.mg.at_node["water__unit_flux_in"])
    with pytest.deprecated_call():
        FlowRouter(dans_grid1.mg, dans_grid1.infile)
    assert_array_equal(np.full(25, 2.0), dans_grid1.mg.at_node["water__unit_flux_in"])


def test_accumulate_D8(dans_grid1):
    """Test accumulation works for D8 in a simple scenario."""
    with pytest.deprecated_call():
        fr = FlowRouter(dans_grid1.mg)
    fr.run_one_step()
    assert_array_equal(dans_grid1.A_target, dans_grid1.mg.at_node["drainage_area"])
    assert_array_equal(
        dans_grid1.frcvr_target, dans_grid1.mg.at_node["flow__receiver_node"]
    )
    assert_array_equal(
        dans_grid1.upids_target, dans_grid1.mg.at_node["flow__upstream_node_order"]
    )
    assert_array_equal(
        dans_grid1.links2rcvr_target,
        dans_grid1.mg.at_node["flow__link_to_receiver_node"],
    )
    assert_array_equal(
        dans_grid1.A_target, dans_grid1.mg.at_node["surface_water__discharge"]
    )
    assert_array_equal(
        dans_grid1.steepest_target, dans_grid1.mg.at_node["topographic__steepest_slope"]
    )


def test_variable_Qin(dans_grid1):
    """Test variable Qin field."""
    Qin_local = np.zeros(25, dtype=float)
    Qin_local[13] = 2.0
    dans_grid1.mg.add_field("node", "water__unit_flux_in", Qin_local, units="m**3/s")
    with pytest.deprecated_call():
        fr = FlowRouter(dans_grid1.mg)
    fr.run_one_step()
    Qout_local = np.zeros_like(Qin_local)
    Qout_local[10:14] = 200.0
    assert_array_equal(Qout_local, dans_grid1.mg.at_node["surface_water__discharge"])
    assert_array_equal(dans_grid1.A_target, dans_grid1.mg.at_node["drainage_area"])
    # note that A DOES NOT CHANGE when messing with Q_in


def test_irreg_topo(dans_grid2):
    """Test D8 routing on a toy irregular topo."""
    with pytest.deprecated_call():
        fr = FlowRouter(dans_grid2.mg)
    fr.run_one_step()
    assert_array_equal(dans_grid2.A_target_D8, dans_grid2.mg.at_node["drainage_area"])
    assert_array_equal(
        dans_grid2.frcvr_target_D8, dans_grid2.mg.at_node["flow__receiver_node"]
    )
    assert_array_equal(
        dans_grid2.upids_target_D8, dans_grid2.mg.at_node["flow__upstream_node_order"]
    )
    assert_array_equal(
        dans_grid2.links2rcvr_target_D8,
        dans_grid2.mg.at_node["flow__link_to_receiver_node"],
    )
    assert dans_grid2.steepest_target_D8 == pytest.approx(
        dans_grid2.mg.at_node["topographic__steepest_slope"]
    )


def test_irreg_topo_old(dans_grid2):
    """
    Test D4 routing on a toy irregular topo, old style, where 'method' is
    passed to the run method, not the init.
    """
    with pytest.deprecated_call():
        fr = FlowRouter(dans_grid2.mg)
    with pytest.deprecated_call():
        fr.route_flow(method="D4")
    assert_array_equal(dans_grid2.A_target_D4, dans_grid2.mg.at_node["drainage_area"])
    assert_array_equal(
        dans_grid2.frcvr_target_D4, dans_grid2.mg.at_node["flow__receiver_node"]
    )
    assert_array_equal(
        dans_grid2.upids_target_D4, dans_grid2.mg.at_node["flow__upstream_node_order"]
    )
    assert_array_equal(
        dans_grid2.links2rcvr_target_D4,
        dans_grid2.mg.at_node["flow__link_to_receiver_node"],
    )
    assert dans_grid2.steepest_target_D4 == pytest.approx(
        dans_grid2.mg.at_node["topographic__steepest_slope"]
    )


def test_irreg_topo_new(dans_grid2):
    """Test D4 routing on a toy irregular topo. 'method' passed to init."""
    with pytest.deprecated_call():
        fr = FlowRouter(dans_grid2.mg, method="D4")
    fr.run_one_step()
    assert_array_equal(dans_grid2.A_target_D4, dans_grid2.mg.at_node["drainage_area"])
    assert_array_equal(
        dans_grid2.frcvr_target_D4, dans_grid2.mg.at_node["flow__receiver_node"]
    )
    assert_array_equal(
        dans_grid2.upids_target_D4, dans_grid2.mg.at_node["flow__upstream_node_order"]
    )
    assert_array_equal(
        dans_grid2.links2rcvr_target_D4,
        dans_grid2.mg.at_node["flow__link_to_receiver_node"],
    )
    assert dans_grid2.steepest_target_D4 == pytest.approx(
        dans_grid2.mg.at_node["topographic__steepest_slope"]
    )


def test_internal_closed(internal_closed):
    """Test closed nodes in the core of the grid."""
    with pytest.deprecated_call():
        fr = FlowRouter(internal_closed.mg)
    fr.run_one_step()
    assert internal_closed.A_target == pytest.approx(
        internal_closed.mg.at_node["drainage_area"]
    )
    assert_array_equal(
        internal_closed.frcvr_target, internal_closed.mg.at_node["flow__receiver_node"]
    )
    assert_array_equal(
        internal_closed.links2rcvr_target,
        internal_closed.mg.at_node["flow__link_to_receiver_node"],
    )
    assert internal_closed.A_target == pytest.approx(
        internal_closed.mg.at_node["surface_water__discharge"]
    )
    assert internal_closed.steepest_target == pytest.approx(
        internal_closed.mg.at_node["topographic__steepest_slope"]
    )


def test_voronoi():
    """Test routing on a (radial) voronoi."""
    vmg = RadialModelGrid(2, dr=2.0)
    z = np.full(20, 10.0, dtype=float)
    # vmg.status_at_node[8:] = CLOSED_BOUNDARY
    all_bounds_but_one = np.array((0, 1, 2, 3, 4, 7, 11, 15, 16, 17, 18, 19))
    vmg.status_at_node[all_bounds_but_one] = CLOSED_BOUNDARY
    # z[7] = 0.  # outlet
    z[12] = 0.0  # outlet
    # inner_elevs = (3., 1., 4., 5., 6., 7., 8.)
    inner_elevs = (8.0, 7.0, 3.0, 1.0, 6.0, 4.0, 5.0)
    # z[:7] = np.array(inner_elevs)
    z[vmg.core_nodes] = np.array(inner_elevs)
    vmg.add_field("node", "topographic__elevation", z, units="-")
    with pytest.deprecated_call():
        fr = FlowRouter(vmg)

    #    nodes_contributing = [np.array([0, 3, 4, 5]),
    #                          np.array([0, 1, 2, 3, 4, 5, 6]),
    #                          np.array([2, ]),
    #                          np.array([3, ]),
    #                          np.array([4, ]),
    #                          np.array([5, ]),
    #                          np.array([6, ])]

    # The follow list contains arrays with the IDs of cells contributing flow
    # to nodes 5, 6, 8, 9, 10, 13, and 14, respectively (which correspond to
    # cells 0-6)
    cells_contributing = [
        np.array([0]),
        np.array([1]),
        np.array([1, 2, 4, 6]),
        np.array([0, 1, 2, 3, 4, 5, 6]),
        np.array([4]),
        np.array([5]),
        np.array([6]),
    ]

    A_target_core = np.zeros(vmg.number_of_core_nodes)
    for i in range(7):
        A_target_core[i] = vmg.area_of_cell[cells_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell.sum()
    fr.run_one_step()
    assert vmg.at_node["drainage_area"][vmg.core_nodes] == pytest.approx(A_target_core)
    assert vmg.at_node["drainage_area"][12] == pytest.approx(A_target_outlet)


def test_voronoi_closedinternal():
    """Test routing on a (radial) voronoi, but with a closed interior node."""
    vmg = RadialModelGrid(2, dr=2.0)
    z = np.full(20, 10.0, dtype=float)
    all_bounds_but_one = np.array((0, 1, 2, 3, 4, 7, 11, 15, 16, 17, 18, 19))
    vmg.status_at_node[all_bounds_but_one] = CLOSED_BOUNDARY
    vmg.status_at_node[8] = CLOSED_BOUNDARY  # new internal closed
    z[12] = 0.0  # outlet
    inner_elevs = (8.0, 7.0, 1.0, 6.0, 4.0, 5.0)
    z[vmg.core_nodes] = np.array(inner_elevs)
    vmg.add_field("node", "topographic__elevation", z, units="-")
    with pytest.deprecated_call():
        fr = FlowRouter(vmg)

    cells_contributing = [
        np.array([0]),
        np.array([1]),
        np.array([0, 1, 3, 4, 5, 6]),
        np.array([1, 4]),
        np.array([1, 4, 5, 6]),
        np.array([1, 4, 6]),
    ]

    A_target_internal = np.zeros(vmg.number_of_core_nodes, dtype=float)
    for i in range(6):
        A_target_internal[i] = vmg.area_of_cell[cells_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell[vmg.cell_at_node[vmg.core_nodes]].sum()
    fr.run_one_step()

    assert vmg.at_node["drainage_area"][vmg.core_nodes] == pytest.approx(
        A_target_internal
    )
    assert vmg.at_node["drainage_area"][12] == pytest.approx(A_target_outlet)
