"""Test the flow accumulator component.

@author: benjaminCampforts
"""

import numpy as np
import pytest
from numpy import testing
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import FieldError
from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import PriorityFloodFlowRouter
from landlab.components.priority_flood_flow_router.cfuncs import _D8_FlowAcc
from landlab.components.priority_flood_flow_router.cfuncs import _D8_flowDir
from landlab.grid.nodestatus import NodeStatus

try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    pytestmark = pytest.mark.skip(reason="richdem is not installed")


def test_check_fields():
    """Check to make sure the right fields have been created."""
    # %%
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    z = mg.add_field("topographic__elevation", mg.node_x**2 + mg.node_y**2, at="node")

    PriorityFloodFlowRouter(mg)
    assert_array_equal(z, mg.at_node["topographic__elevation"])
    assert_array_equal(mg.at_node["drainage_area"], 0.0)
    assert_array_equal(mg.at_node["water__unit_flux_in"], 1.0)

    PriorityFloodFlowRouter(mg, runoff_rate=2.0)
    assert_array_equal(mg.at_node["water__unit_flux_in"], 2.0)


def input_values():
    """
    PriorityFloodFlowRouter should throw an error when wrong input values are provided
    """

    # %% Default configuration
    grid = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    grid.add_field("topographic__elevation", grid.node_x + grid.node_y, at="node")

    with pytest.raises(ValueError):
        PriorityFloodFlowRouter(grid, suppress_out=True, flow_metric="Oops")

    with pytest.raises(ValueError):
        PriorityFloodFlowRouter(grid, suppress_out=True, depression_handler="Oops")

    assert "water__unit_flux_in" not in grid.at_node
    PriorityFloodFlowRouter(grid)
    assert_array_almost_equal(grid.at_node["water__unit_flux_in"], 1.0)

    grid.at_node.pop("water__unit_flux_in")
    assert "water__unit_flux_in" not in grid.at_node
    PriorityFloodFlowRouter(grid, runoff_rate=999)
    assert_array_almost_equal(grid.at_node["water__unit_flux_in"], 999)

    grid.at_node["water__unit_flux_in"].fill(10.0)
    assert "water__unit_flux_in" in grid.at_node
    PriorityFloodFlowRouter(grid)
    assert_array_almost_equal(grid.at_node["water__unit_flux_in"], 10.0)

    assert "water__unit_flux_in" in grid.at_node
    PriorityFloodFlowRouter(grid, runoff_rate=5.0)
    assert_array_almost_equal(grid.at_node["water__unit_flux_in"], 5.0)


# %%
def test_seperate_hillflow_update():
    # %%
    """
    Hillflow fields cannot be updated if not initialized
    """

    # %% Default configuration
    mg1 = RasterModelGrid((6, 6), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x + mg1.node_y, at="node")
    fa1 = PriorityFloodFlowRouter(mg1, suppress_out=True)
    fa1.run_one_step()
    with pytest.raises(ValueError):
        fa1.update_hill_fdfa()

    # %%


def test_fields1():
    # %%
    """
    Check to make sure the right fields have been created.

    Check that the sizes are also correct.
    """

    # %% Default configuration
    mg1 = RasterModelGrid((6, 6), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x + mg1.node_y, at="node")
    fa1 = PriorityFloodFlowRouter(mg1, suppress_out=True)
    fa1.run_one_step()
    assert sorted(mg1.at_node) == [
        "depression_free_elevation",
        "drainage_area",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "squared_length_adjacent",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]


def test_fields2():
    # %% No flow accumulation
    mg2 = RasterModelGrid((7, 7), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x + mg2.node_y, at="node")
    fa2 = PriorityFloodFlowRouter(mg2, suppress_out=True, accumulate_flow=False)
    fa2.run_one_step()

    assert sorted(mg2.at_node) == [
        "depression_free_elevation",
        "flood_status_code",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__upstream_node_order",
        "squared_length_adjacent",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]
    # %%


def test_fields3():
    # %% Second FD (no FA is default)
    mg3 = RasterModelGrid((8, 8), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x + mg3.node_y, at="node")
    fa3 = PriorityFloodFlowRouter(mg3, separate_hill_flow=True, suppress_out=True)
    fa3.run_one_step()
    assert sorted(mg3.at_node) == [
        "depression_free_elevation",
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
        "squared_length_adjacent",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]
    # %%


def test_fields4():
    """
    Second FD (with FA )
    """
    mg4 = RasterModelGrid((9, 9), xy_spacing=(1, 1))
    mg4.add_field("topographic__elevation", mg4.node_x + mg4.node_y, at="node")
    fa4 = PriorityFloodFlowRouter(
        mg4, separate_hill_flow=True, accumulate_flow_hill=True, suppress_out=True
    )
    fa4.run_one_step()
    assert sorted(mg4.at_node) == [
        "depression_free_elevation",
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
        "squared_length_adjacent",
        "surface_water__discharge",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]


def test_fields5():
    # %% Verify multiple flow
    mg5 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg5.add_field("topographic__elevation", mg5.node_x + mg5.node_y, at="node")
    _ = PriorityFloodFlowRouter(mg5, suppress_out=True, flow_metric="Quinn")

    assert mg5.at_node["topographic__steepest_slope"].shape[1] == 8
    assert mg5.at_node["flow__receiver_node"].shape[1] == 8
    assert mg5.at_node["flow__receiver_proportions"].shape[1] == 8
    assert mg5.at_node["flow__link_to_receiver_node"].shape[1] == 8


def test_accumulated_area_closes():
    """Check that accumulated area is area of core nodes."""
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(mg, suppress_out=True)
    fa.run_one_step()

    drainage_area = mg.at_node["drainage_area"]
    drained_area = np.sum(drainage_area[mg.boundary_nodes])
    core_area = np.sum(mg.cell_area_at_node[mg.core_nodes])
    assert drained_area == core_area

    # %% multiple flow
    """Check that accumulated area is area of core nodes."""
    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x + mg3.node_y, at="node")
    fa3 = PriorityFloodFlowRouter(
        mg3,
        separate_hill_flow=True,
        accumulate_flow_hill=True,
        hill_flow_metric="Quinn",
        suppress_out=True,
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
        PriorityFloodFlowRouter(mg, flow_director="D2")


# %%
def test_D8_metric():
    # %%
    """Test D8 routing functionality"""

    topographic__elevation = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 21.0, 10.0, 0.0],
            [0.0, 31.0, 30.8, 0.0],
            [0.0, 32.0, 30.9, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )

    mg2 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", topographic__elevation, at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    fa2 = PriorityFloodFlowRouter(mg2, flow_metric="D8", suppress_out=True)
    fa2.run_one_step()
    pf_r = mg2.at_node["flow__receiver_node"]
    reciever = np.array(
        [
            [0, 1, 2, 3],
            [4, 1, 2, 7],
            [8, 6, 6, 11],
            [12, 14, 10, 15],
            [16, 17, 18, 19],
        ]
    )
    assert_array_equal(reciever.flatten(), pf_r)

    # Compare with default Landlab flow accumulator
    mg1 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    from landlab.components.flow_accum import FlowAccumulator

    mg1.add_field("topographic__elevation", topographic__elevation, at="node")
    mg1.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa1 = FlowAccumulator(mg1, flow_director="D8")
    fa1.run_one_step()
    defaultLandl_r = mg1.at_node["flow__receiver_node"]
    assert_array_equal(defaultLandl_r, pf_r)


# %%
def test_various_surface_values():
    # %%

    topographic__elevation = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 30.9, 21.0, 0.0],
            [0.0, 30.8, 31.8, 0.0],
            [0.0, 10.0, 11.1, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )

    topography = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 121.0, 110.0, 0.0],
            [0.0, 131.0, 130.8, 0.0],
            [0.0, 132.0, 130.9, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )

    mg2 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", topographic__elevation, at="node")
    mg2.add_field("topography", topography, at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa2 = PriorityFloodFlowRouter(
        mg2, flow_metric="D8", suppress_out=True, surface="topography"
    )
    fa2.run_one_step()
    pf_r = mg2.at_node["flow__receiver_node"]
    reciever = np.array(
        [
            [0, 1, 2, 3],
            [4, 1, 2, 7],
            [8, 6, 6, 11],
            [12, 14, 10, 15],
            [16, 17, 18, 19],
        ]
    )
    assert_array_equal(reciever.flatten(), pf_r)

    fa3 = PriorityFloodFlowRouter(
        mg2, flow_metric="D8", suppress_out=True, surface=topography
    )
    fa3.run_one_step()
    pf_r = mg2.at_node["flow__receiver_node"]
    assert_array_equal(reciever.flatten(), pf_r)

    topo_fd = fa3.surface_values
    assert_array_equal(topo_fd, np.reshape(topography, 4 * 5))


# %%
def test_flow_accumulator_properties():
    # %%

    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(mg, flow_metric="D8", suppress_out=True)
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            [3.0, 2.0, 1.0, 0.0, 0.0],
            [2.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 2.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    assert_array_equal(fa.node_water_discharge, node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())


def test_flow_accumulator_properties_D4():
    # %% set update_hill_flow_instantaneous to false and update hill flow properties explicitly

    mg = RasterModelGrid((4, 4), xy_spacing=(1, 1))
    _ = mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(
        mg,
        flow_metric="D4",
        suppress_out=True,
        separate_hill_flow=True,
        hill_flow_metric="D4",
        accumulate_flow_hill=True,
        update_hill_flow_instantaneous=False,
    )
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()
    fa.update_hill_fdfa()

    node_drainage_area = np.array(
        [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 1.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )

    assert_array_equal(fa.node_water_discharge, node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())

    # Hill flow accumulation should be similar
    assert_array_equal(
        mg.at_node["hill_surface_water__discharge"], node_drainage_area.flatten()
    )
    assert_array_equal(mg.at_node["hill_drainage_area"], node_drainage_area.flatten())


# %%
def test_flow_accumulator_properties_D4_varying_runoff():
    # %%

    mg = RasterModelGrid((4, 4), xy_spacing=(1, 1))
    _ = mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    runoff_rate = mg.node_x * 0 + 2
    fa = PriorityFloodFlowRouter(
        mg,
        flow_metric="D4",
        suppress_out=True,
        separate_hill_flow=True,
        hill_flow_metric="D4",
        runoff_rate=runoff_rate,
        accumulate_flow_hill=True,
    )
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 1.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )

    assert_array_equal(fa.node_water_discharge, 2 * node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())

    # Hill flow accumulation should be similar
    assert_array_equal(
        mg.at_node["hill_surface_water__discharge"], 2 * node_drainage_area.flatten()
    )
    assert_array_equal(mg.at_node["hill_drainage_area"], node_drainage_area.flatten())


def test_flow_accumulator_properties_D4_varying_water__unit_flux_in():
    # %%

    mg = RasterModelGrid((4, 4), xy_spacing=(1, 1))
    _ = mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    mg.add_field("water__unit_flux_in", mg.node_x * 0 + 2, at="node")
    fa = PriorityFloodFlowRouter(
        mg,
        flow_metric="D4",
        suppress_out=True,
        separate_hill_flow=True,
        hill_flow_metric="D4",
        accumulate_flow_hill=True,
    )
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 1.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )

    assert_array_equal(fa.node_water_discharge, 2 * node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())

    # Hill flow accumulation should be similar
    assert_array_equal(
        mg.at_node["hill_surface_water__discharge"], 2 * node_drainage_area.flatten()
    )
    assert_array_equal(mg.at_node["hill_drainage_area"], node_drainage_area.flatten())


# %%
def test_flow_accumulator_varying_Runoff_rate():
    """
    There are two options to provide a runoff value
    The first one is to provide a weight array as a runoff_value input argument
    """

    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    runoff_rate = mg.node_x * 0 + 2
    fa = PriorityFloodFlowRouter(
        mg,
        flow_metric="D8",
        suppress_out=True,
        runoff_rate=runoff_rate,
        separate_hill_flow=True,
        hill_flow_metric="D8",
        accumulate_flow_hill=True,
    )

    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            [3.0, 2.0, 1.0, 0.0, 0.0],
            [2.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 2.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    assert_array_equal(fa.node_water_discharge, 2 * node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())

    # Hill flow accumulation should be similar
    assert_array_equal(
        mg.at_node["hill_surface_water__discharge"], 2 * node_drainage_area.flatten()
    )
    assert_array_equal(mg.at_node["hill_drainage_area"], node_drainage_area.flatten())


# %%
def test_flow_accumulator_varying_water__unit_flux_in():
    """
    There are two options to provide a runoff value
    The second is to provide a weight array as a water__unit_flux_in input field
    """

    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    mg.add_field("water__unit_flux_in", mg.node_x * 0 + 2, at="node")
    fa = PriorityFloodFlowRouter(
        mg,
        flow_metric="D8",
        suppress_out=True,
        separate_hill_flow=True,
        hill_flow_metric="D8",
        accumulate_flow_hill=True,
    )

    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            [3.0, 2.0, 1.0, 0.0, 0.0],
            [2.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 2.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    assert_array_equal(fa.node_water_discharge, 2 * node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())

    # Hill flow accumulation should be similar
    assert_array_equal(
        mg.at_node["hill_surface_water__discharge"], 2 * node_drainage_area.flatten()
    )
    assert_array_equal(mg.at_node["hill_drainage_area"], node_drainage_area.flatten())


# %%
def test_flow_accumulator_properties_breach():
    # %%

    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x * 2 + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(
        mg, flow_metric="D8", suppress_out=True, depression_handler="breach"
    )
    # from landlab.components.flow_accum import FlowAccumulator
    # fa = FlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            [3.0, 2.0, 1.0, 0.0, 0.0],
            [2.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 2.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    assert_array_equal(fa.node_water_discharge, node_drainage_area.flatten())
    assert_array_equal(fa.node_drainage_area, node_drainage_area.flatten())


# %%
def test_bad_metric_name():
    # %%
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        PriorityFloodFlowRouter(mg, flow_metric="spam", suppress_out=True)


# %%
def test_bad_hill_metric_name():
    # %%
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        PriorityFloodFlowRouter(
            mg, hill_flow_metric="Landlab_is_cool", suppress_out=True
        )


# %%
def test_bad_depression_handler():
    # %%
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        PriorityFloodFlowRouter(
            mg, depression_handler="depression_filler", suppress_out=True
        )


# %%
def test_hex_mfd():
    # %%
    mg = HexModelGrid((5, 3))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(FieldError):
        PriorityFloodFlowRouter(mg, flow_metric="MFD", suppress_out=True)
    # %%


def test_sum_prop_is_one():
    # %%
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    fa.run_one_step()

    # Single flow
    props_Pf = mg.at_node["flow__receiver_proportions"]
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow proportions is not equal to one",
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
        err_msg="Sum of flow proportions is not equal to one",
        verbose=True,
    )

    # %% multiple flow with D8 over hills
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(
        mg, separate_hill_flow=True, hill_flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Multiple flow
    props_Pf = mg.at_node["hill_flow__receiver_proportions"]
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow proportions is not equal to one",
        verbose=True,
    )

    # %% multiple flow with D8 over rivers
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = PriorityFloodFlowRouter(mg, flow_metric="Quinn", suppress_out=True)
    fa.run_one_step()

    # Multiple flow
    props_Pf = mg.at_node["flow__receiver_proportions"]
    props_Pf[props_Pf == -1] = 0
    props_Pf = np.sum(props_Pf, axis=1)
    testing.assert_array_almost_equal(
        props_Pf,
        np.ones_like(props_Pf),
        decimal=5,
        err_msg="Sum of flow proportions is not equal to one",
        verbose=True,
    )
    # %%


def test_cython_functions():
    # %% test d8 flow dir with open boundaries
    mg = RasterModelGrid((4, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[5:7] = [1, 3]
    z[9:11] = [3, 4]
    from landlab.plot import imshow_grid

    imshow_grid(mg, "topographic__elevation")

    activeCells = np.array(mg.status_at_node != NodeStatus.CLOSED + 0, dtype=int)
    receivers = np.array(mg.status_at_node, dtype=int)
    distance_receiver = np.zeros((receivers.shape), dtype=float)
    cores = mg.core_nodes
    activeCores = cores[activeCells[cores] == 1]
    # Make boundaries to save time with conditionals in c loops
    receivers[np.nonzero(mg.status_at_node)] = -1
    steepest_slope = np.zeros((receivers.shape), dtype=float)
    el_dep_free = z
    el_ori = mg.at_node["topographic__elevation"]
    dist = (
        np.array([1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2)]) * mg.dx
    )

    ngb = np.zeros((8,), dtype=int)
    el_d = np.zeros((8,), dtype=float)

    # Links
    adj_link = np.array(mg.d8s_at_node, dtype=int)
    recvr_link = np.zeros((receivers.shape), dtype=int) - 1

    _D8_flowDir(
        receivers,
        distance_receiver,
        steepest_slope,
        np.array(el_dep_free),
        el_ori,
        dist,
        ngb,
        activeCores,
        activeCells,
        el_d,
        mg.number_of_node_columns,
        mg.dx,
        adj_link,
        recvr_link,
    )

    # We know where the water will flow using the D8 steepest descent algo
    # Also consider that under equal slopes, the flow will follow Landlab's
    # rotational ordering going first to cardial, then to diagonal cells
    known_rec = np.array([-1, -1, -1, -1, -1, 4, 7, -1, -1, 13, 11, -1, -1, -1, -1, -1])

    testing.assert_array_equal(
        known_rec,
        receivers,
        err_msg="Error with D8 flow routing calculations",
        verbose=True,
    )

    # %% test flow acc with open boundaries
    node_cell_area = np.array(mg.cell_area_at_node)
    node_cell_area[mg.closed_boundary_nodes] = 0.0
    dis = np.full(mg.number_of_nodes, node_cell_area)

    da = np.array(node_cell_area)
    sort = np.argsort(el_dep_free)
    stack_flip = np.flip(sort)
    # Filter out donors giving to receivers being -1
    stack_flip = np.array(stack_flip[receivers[stack_flip] != -1], dtype=int)

    _D8_FlowAcc(da, dis, stack_flip, receivers)

    # We know how much water will accumualte given the calculated flow dir (see before)
    known_FA = np.array(
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0]
    )

    testing.assert_array_equal(
        known_FA,
        da,
        err_msg="Error with D8_FlowAcc calculations",
        verbose=True,
    )

    # %% test d8 flow dir with closed boundaries
    mg = RasterModelGrid((4, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    # Close all model boundary edges
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # Set lower-left (southwest) corner as an open boundary
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )

    z[5:7] = [1, 3]
    z[9:11] = [3, 4]

    activeCells = np.array(mg.status_at_node != NodeStatus.CLOSED + 0, dtype=int)
    receivers = np.array(mg.status_at_node, dtype=int)
    distance_receiver = np.zeros((receivers.shape), dtype=float)
    cores = mg.core_nodes
    activeCores = cores[activeCells[cores] == 1]
    # Make boundaries to save time with conditionals in c loops
    receivers[np.nonzero(mg.status_at_node)] = -1
    steepest_slope = np.zeros((receivers.shape), dtype=float)
    el_dep_free = z
    el_ori = mg.at_node["topographic__elevation"]
    dist = (
        np.array([1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2)]) * mg.dx
    )

    ngb = np.zeros((8,), dtype=int)
    el_d = np.zeros((8,), dtype=float)

    # Links
    adj_link = np.array(mg.d8s_at_node, dtype=int)
    recvr_link = np.zeros((receivers.shape), dtype=int) - 1

    _D8_flowDir(
        receivers,
        distance_receiver,
        steepest_slope,
        np.array(el_dep_free),
        el_ori,
        dist,
        ngb,
        activeCores,
        activeCells,
        el_d,
        mg.number_of_node_columns,
        mg.dx,
        adj_link,
        recvr_link,
    )

    # We know where the water will flow using the D8 steepest descent algo
    # Also consider that under equal slopes, the flow will follow Landlab's
    # rotational ordering going first to cardial, then to diagonal cells
    known_rec = np.array([-1, -1, -1, -1, -1, 0, 5, -1, -1, 5, 5, -1, -1, -1, -1, -1])

    testing.assert_array_equal(
        known_rec,
        receivers,
        err_msg="Error with D8 flow routing calculations",
        verbose=True,
    )

    # %% test flow acc with closed boundaries
    node_cell_area = np.array(mg.cell_area_at_node)
    node_cell_area[mg.closed_boundary_nodes] = 0.0
    dis = np.full(mg.number_of_nodes, node_cell_area)

    da = np.array(node_cell_area)
    sort = np.argsort(el_dep_free)
    stack_flip = np.flip(sort)
    # Filter out donors giving to receivers being -1
    stack_flip = np.array(stack_flip[receivers[stack_flip] != -1], dtype=int)

    _D8_FlowAcc(da, dis, stack_flip, receivers)

    # We know how much water will accumualte given the calculated flow dir (see before)
    known_FA = np.array(
        [4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )

    testing.assert_array_equal(
        known_FA,
        da,
        err_msg="Error with D8_FlowAcc calculations",
        verbose=True,
    )


# %% test flooded nodes
def test_flooded_nodes():
    mg = RasterModelGrid((6, 6), xy_spacing=(1, 1))
    z = mg.add_ones("topographic__elevation", at="node")
    z[mg.boundary_nodes] = 0
    z[14:16] = 0.5
    z[20:22] = 0.75
    pf = PriorityFloodFlowRouter(mg)
    pf.run_one_step()
    # FLOODED nodes have code 3
    testing.assert_array_equal(12, np.sum(mg.at_node["flood_status_code"]))


# %% test boundary conditions
def test_boundary_conditions():
    mg = RasterModelGrid((10, 10), 1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.at_node["topographic__elevation"][mg.core_nodes] = np.random.rand(
        len(mg.core_nodes)
    )
    # set boundary conditions - right, top, left, bottom
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    # All of the interior nodes flow to outlet node number 5
    mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    flow = PriorityFloodFlowRouter(mg, flow_metric="D8", update_flow_depressions=True)
    flow.run_one_step()
    # on a 10 by 10 grid, there are 8 by 8 (64) cells draining to the outlet node (5)
    # if boundary conditions are properly set
    testing.assert_array_equal(64, mg.at_node["drainage_area"][5])
