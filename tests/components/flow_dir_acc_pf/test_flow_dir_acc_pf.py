"""Test the flow accumulator component.

@author: benjaminCampforts
"""
import numpy as np
import pytest
from numpy import testing
from numpy.testing import assert_array_equal

from landlab import FieldError, HexModelGrid, RasterModelGrid
from landlab.components import FlowDirAccPf


def test_check_fields():
    """Check to make sure the right fields have been created."""
    # %%
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    z = mg.add_field(
        "topographic__elevation", mg.node_x ** 2 + mg.node_y ** 2, at="node"
    )

    FlowDirAccPf(mg)
    assert_array_equal(z, mg.at_node["topographic__elevation"])
    assert_array_equal(np.zeros(100), mg.at_node["drainage_area"])
    assert_array_equal(np.ones(100), mg.at_node["water__unit_flux_in"])

    FlowDirAccPf(mg, runoff_rate=2.0)
    assert_array_equal(np.full(100, 2.0), mg.at_node["water__unit_flux_in"])


# %%
def test_fields():
    # %%
    """Check to make sure the right fields have been created.

    Check that the sizes are also correct.
    """

    # %% Default configuration
    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x + mg1.node_y, at="node")
    fa1 = FlowDirAccPf(mg1)
    fa1.run_one_step()

    assert sorted(list(mg1.at_node.keys())) == [
        "SQUARED_length_adjacent",
        "deprFree_elevation",
        "drainage_area",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]
    # %% No flow accumulation
    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x + mg2.node_y, at="node")
    fa2 = FlowDirAccPf(mg2, accumulate_flow=True)
    fa2.run_one_step()
    assert sorted(list(mg2.at_node.keys())) == [
        "SQUARED_length_adjacent",
        "deprFree_elevation",
        "drainage_area",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]

    # %% Second FD (no FA is default)
    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x + mg2.node_y, at="node")
    fa2 = FlowDirAccPf(mg2, separate_hill_flow=True)
    fa2.run_one_step()
    assert sorted(list(mg2.at_node.keys())) == [
        "SQUARED_length_adjacent",
        "deprFree_elevation",
        "drainage_area",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "hill_flow__receiver_node",
        "hill_flow__receiver_proportions",
        "hill_flow__upstream_node_order",
        "hill_topographic__steepest_slope",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]
    # %% Second FD (with FA )
    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x + mg2.node_y, at="node")
    fa2 = FlowDirAccPf(mg2, separate_hill_flow=True, accumulate_flow_hill=True)
    fa2.run_one_step()
    assert sorted(list(mg2.at_node.keys())) == [
        "SQUARED_length_adjacent",
        "deprFree_elevation",
        "drainage_area",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "hill_drainage_area",
        "hill_flow__receiver_node",
        "hill_flow__receiver_proportions",
        "hill_flow__upstream_node_order",
        "hill_surface_water__discharge",
        "hill_topographic__steepest_slope",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]

    # %%


def test_accumulated_area_closes():
    """Check that accumulated area is area of core nodes."""
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = FlowDirAccPf(mg)
    fa.run_one_step()

    drainage_area = mg.at_node["drainage_area"]
    drained_area = np.sum(drainage_area[mg.boundary_nodes])
    core_area = np.sum(mg.cell_area_at_node[mg.core_nodes])
    assert drained_area == core_area

    # %% multiple flow
    """Check that accumulated area is area of core nodes."""
    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x + mg3.node_y, at="node")
    fa3 = FlowDirAccPf(
        mg3,
        separate_hill_flow=True,
        accumulate_flow_hill=True,
        hill_flow_metric="Quinn",
    )
    fa3.run_one_step()

    drainage_area_hill = mg3.at_node["hill_drainage_area"]
    drainage_area_hill = np.sum(drainage_area_hill[mg3.boundary_nodes])
    core_area = np.sum(mg3.cell_area_at_node[mg3.core_nodes])
    assert drainage_area_hill == core_area


# %%


def test_specifying_routing_method_wrong():
    # %%
    """Test specifying incorrect method for routing compatability"""
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    with pytest.raises(TypeError):
        FlowDirAccPf(mg, flow_director="D2")


# %%
def test_D8_metric():
    # %%
    """Test D8 routing functionality"""

    topographic__elevation = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            21.0,
            10.0,
            0.0,
            0.0,
            31.0,
            30.8,
            0.0,
            0.0,
            32.0,
            30.9,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    mg2 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", topographic__elevation, at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    fa2 = FlowDirAccPf(mg2, flow_metric="D8")
    fa2.run_one_step()
    pf_r = mg2.at_node["flow__receiver_node"]
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 6, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )
    assert_array_equal(reciever, pf_r)

    # Uncomment to compare with default Landlab flow accumulator
    mg1 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    from landlab.components.flow_accum import FlowAccumulator

    mg1.add_field("topographic__elevation", topographic__elevation, at="node")
    mg1.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa1 = FlowAccumulator(mg1, flow_director="D8")
    fa1.run_one_step()
    defaultLandl_r = mg1.at_node["flow__receiver_node"]
    assert_array_equal(defaultLandl_r, pf_r)


# %%
def test_flow_accumulator_properties():
    # %%

    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    fa = FlowDirAccPf(mg, flow_metric="D8")
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            3.0,
            2.0,
            1.0,
            0.0,
            0.0,
            2.0,
            3.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            2.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(fa.node_water_discharge, node_drainage_area)
    assert_array_equal(fa.node_drainage_area, node_drainage_area)


# %%
def test_bad_metric_name():
    # %%
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(Exception):
        FlowDirAccPf(mg, flow_metric="spam")


# %%
def test_hex_mfd():
    # %%
    mg = HexModelGrid((5, 3))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(FieldError):
        FlowDirAccPf(mg, flow_metric="MFD")
    # %%


def test_sum_prop_is_one():
    # %%
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = FlowDirAccPf(mg, separate_hill_flow=True)
    fa.run_one_step()

    # Single flow
    props_Pf = mg.at_node["flow__receiver_proportions"]
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow propoertions is not equal to one",
        verbose=True,
    )

    # Multiple flow
    props_Pf = mg.at_node["hill_flow__receiver_proportions"]
    props_Pf[props_Pf == -1] = 0
    props_Pf = props_Pf.sum(axis=1)
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow propoertions is not equal to one",
        verbose=True,
    )

    # %% multiple flow with D8
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = FlowDirAccPf(mg, separate_hill_flow=True, hill_flow_metric="D8")
    fa.run_one_step()

    # Multiple flow
    props_Pf = mg.at_node["hill_flow__receiver_proportions"]
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow propoertions is not equal to one",
        verbose=True,
    )
