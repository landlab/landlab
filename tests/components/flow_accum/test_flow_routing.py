import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RadialModelGrid
from landlab.components import FlowAccumulator

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_check_fields(dans_grid1):
    """Check to make sure the right fields have been created."""
    FlowAccumulator(dans_grid1.mg, flow_director="D8")
    assert_array_equal(dans_grid1.z, dans_grid1.mg.at_node["topographic__elevation"])
    assert_array_equal(np.zeros(25), dans_grid1.mg.at_node["drainage_area"])
    assert_array_equal(np.ones(25), dans_grid1.mg.at_node["water__unit_flux_in"])
    FlowAccumulator(dans_grid1.mg, runoff_rate=2.0)
    assert_array_equal(np.full(25, 2.0), dans_grid1.mg.at_node["water__unit_flux_in"])


def test_check_field_input(dans_grid1):
    """Check we can successfully pass water__discharge_in."""
    dans_grid1.mg.add_field(
        "water__unit_flux_in", np.full(25, 3.0), at="node", units="m**3/s"
    )
    FlowAccumulator(dans_grid1.mg, flow_director="D8")

    assert_array_equal(np.full(25, 3.0), dans_grid1.mg.at_node["water__unit_flux_in"])
    FlowAccumulator(dans_grid1.mg, runoff_rate=2.0)
    assert_array_equal(np.full(25, 2.0), dans_grid1.mg.at_node["water__unit_flux_in"])


def test_accumulate_D8(dans_grid1):
    """Test accumulation works for D8 in a simple scenario."""
    fr = FlowAccumulator(dans_grid1.mg, flow_director="D8")

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
    dans_grid1.mg.add_field("water__unit_flux_in", Qin_local, at="node", units="m**3/s")
    fr = FlowAccumulator(dans_grid1.mg, flow_director="D8")

    fr.run_one_step()
    Qout_local = np.zeros_like(Qin_local)
    Qout_local[10:14] = 200.0
    assert_array_equal(Qout_local, dans_grid1.mg.at_node["surface_water__discharge"])
    assert_array_equal(dans_grid1.A_target, dans_grid1.mg.at_node["drainage_area"])


def test_irreg_topo(dans_grid2):
    """Test D8 routing on a toy irregular topo."""
    fr = FlowAccumulator(dans_grid2.mg, flow_director="D8")

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


def test_irreg_topo_new(dans_grid2):
    """Test D4 routing on a toy irregular topo. 'method' passed to init."""
    fr = FlowAccumulator(dans_grid2.mg, flow_director="D4")
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
    fr = FlowAccumulator(internal_closed.mg, flow_director="D8")

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
    vmg = RadialModelGrid(n_rings=2, nodes_in_first_ring=6)

    outlet_node = 11
    z = np.full(vmg.number_of_nodes, 10.0)

    all_bounds_but_one = [0, 1, 2, 3, 4, 7, 14, 15, 16, 17, 18]
    vmg.status_at_node[all_bounds_but_one] = vmg.BC_NODE_IS_CLOSED
    z[outlet_node] = 0.0  # outlet
    z[vmg.core_nodes] = [7.0, 8.0, 6.0, 3.0, 1.0, 5.0, 4.0]
    vmg.add_field("topographic__elevation", z, at="node", units="-")
    fr = FlowAccumulator(vmg)

    # The follow list contains arrays with the IDs of cells contributing flow
    # to nodes 5, 6, 8, 9, 10, 13, and 14, respectively (which correspond to
    # cells 0-6)
    cells_contributing = [
        np.array([0]),
        np.array([1]),
        np.array([2]),
        np.array([0, 3, 2, 5]),
        np.array([0, 1, 2, 3, 4, 5, 6]),
        np.array([5]),
        np.array([6]),
    ]

    A_target_core = np.zeros(len(vmg.core_nodes))
    for i in range(len(cells_contributing)):
        A_target_core[i] = vmg.area_of_cell[cells_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell.sum()
    fr.run_one_step()

    assert vmg.at_node["drainage_area"][vmg.core_nodes] == pytest.approx(A_target_core)
    assert vmg.at_node["drainage_area"][11] == pytest.approx(A_target_outlet)


def test_voronoi_closedinternal():
    """Test routing on a (radial) voronoi, but with a closed interior node."""
    vmg = RadialModelGrid(n_rings=2, nodes_in_first_ring=6)

    outlet_node = 11
    z = np.full(vmg.number_of_nodes, 10.0)

    all_bounds_but_one = [0, 1, 2, 3, 4, 7, 14, 15, 16, 17, 18]
    vmg.status_at_node[all_bounds_but_one] = vmg.BC_NODE_IS_CLOSED
    vmg.status_at_node[9] = vmg.BC_NODE_IS_CLOSED  # new internal closed
    z[outlet_node] = 0.0
    z[vmg.core_nodes] = [7.0, 8.0, 6.0, 1.0, 5.0, 4.0]
    vmg.add_field("topographic__elevation", z, at="node", units="-")
    fr = FlowAccumulator(vmg)

    cells_contributing = [
        np.array([0]),
        np.array([1]),
        np.array([0, 2]),
        np.array([0, 1, 2, 4, 5, 6]),
        np.array([0, 2, 5]),
        np.array([0, 2, 5, 6]),
    ]

    A_target_internal = np.zeros(len(vmg.core_nodes))
    for i in range(len(cells_contributing)):
        A_target_internal[i] = vmg.area_of_cell[cells_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell[vmg.cell_at_node[vmg.core_nodes]].sum()
    fr.run_one_step()

    assert vmg.at_node["drainage_area"][vmg.core_nodes] == pytest.approx(
        A_target_internal
    )
    assert vmg.at_node["drainage_area"][outlet_node] == pytest.approx(A_target_outlet)
