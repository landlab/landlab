"""Test the flow accumulator component.

@author: krb
"""

# Created on Thurs Nov 12, 2015
import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import LinearDiffuser
from landlab.components.depression_finder.lake_mapper import DepressionFinderAndRouter
from landlab.components.flow_accum import LossyFlowAccumulator
from landlab.components.flow_director import FlowDirectorD8
from landlab.components.flow_director import FlowDirectorDINF
from landlab.components.flow_director import FlowDirectorMFD
from landlab.components.flow_director import FlowDirectorSteepest

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_loss_func_arguments():
    """Check the loss_function only has one argument."""
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x**2 + mg.node_y**2, at="node")

    def funcgood(x):
        return 0.0

    def funcgood2(x, y):
        return 0.0

    def funcgood3(x, y, z):
        return 0.0

    def funcgood4(x, y, z, grid):
        return 0.0

    def funcbad(x, y, z, grid, K):
        return 0.0

    def funcbad2(x):
        return "booooooo"

    def funcbad3(x, y):
        return np.array([0.0])

    def funcbad4(x, y, z):
        return (0, 1)

    LossyFlowAccumulator(mg, loss_function=funcgood)  # no problem
    LossyFlowAccumulator(mg, loss_function=funcgood2)  # no problem
    LossyFlowAccumulator(mg, loss_function=funcgood3)  # no problem
    LossyFlowAccumulator(mg, loss_function=funcgood4)  # no problem
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, loss_function=funcbad)
    with pytest.raises(TypeError):
        LossyFlowAccumulator(mg, loss_function=funcbad2)
    with pytest.raises(TypeError):
        LossyFlowAccumulator(mg, loss_function=funcbad3)
    with pytest.raises(TypeError):
        LossyFlowAccumulator(mg, loss_function=funcbad4)


def test_run_with_2_fn_args():
    mg = RasterModelGrid((3, 5), xy_spacing=(2, 1))
    mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    def mylossfunction(qw, nodeID):
        return 0.5 * qw

    fa = LossyFlowAccumulator(
        mg,
        "topographic__elevation",
        flow_director=FlowDirectorSteepest,
        loss_function=mylossfunction,
    )
    fa.run_one_step()

    assert np.allclose(
        mg.at_node["drainage_area"].reshape(mg.shape),
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [6.0, 6.0, 4.0, 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ),
    )
    assert np.allclose(
        mg.at_node["surface_water__discharge"].reshape(mg.shape),
        np.array(
            [
                [0.00, 0.00, 0.00, 0.00, 0.00],
                [1.75, 3.50, 3.00, 2.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00],
            ]
        ),
    )


def test_run_with_3_fn_args():
    mg = RasterModelGrid((3, 5), xy_spacing=(2, 1))
    mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    def mylossfunction(qw, nodeID, linkID):
        return 0.5 * qw

    fa = LossyFlowAccumulator(
        mg,
        "topographic__elevation",
        flow_director=FlowDirectorSteepest,
        loss_function=mylossfunction,
    )
    fa.run_one_step()

    assert np.allclose(
        mg.at_node["drainage_area"].reshape(mg.shape),
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [6.0, 6.0, 4.0, 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ),
    )
    assert np.allclose(
        mg.at_node["surface_water__discharge"].reshape(mg.shape),
        np.array(
            [
                [0.00, 0.00, 0.00, 0.00, 0.00],
                [1.75, 3.50, 3.00, 2.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00],
            ]
        ),
    )


def test_check_fields():
    """Check to make sure the right fields have been created."""

    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    z = mg.add_field("topographic__elevation", mg.node_x**2 + mg.node_y**2, at="node")

    LossyFlowAccumulator(mg)
    assert_array_equal(z, mg.at_node["topographic__elevation"])
    assert_array_equal(np.zeros(100), mg.at_node["drainage_area"])
    assert_array_equal(np.ones(100), mg.at_node["water__unit_flux_in"])
    assert_array_equal(np.zeros(100), mg.at_node["surface_water__discharge_loss"])

    LossyFlowAccumulator(mg, runoff_rate=2.0)
    assert_array_equal(np.full(100, 2.0), mg.at_node["water__unit_flux_in"])

    # quick test that the component binds correctly to an existing field:
    L = mg.at_node["surface_water__discharge_loss"]
    fa = LossyFlowAccumulator(mg)
    fa.run_one_step()  # this line is padding to make flake8 happy
    L[0] = 1.0
    assert mg.at_node["surface_water__discharge_loss"] is L


def test_director_adding_methods_are_equivalent_Steepest():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg0.add_field("topographic__elevation", mg0.node_x**2 + mg0.node_y**2, at="node")
    fa0 = LossyFlowAccumulator(mg0, flow_director="D4")
    fa0.run_one_step()

    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x**2 + mg1.node_y**2, at="node")
    fa1 = LossyFlowAccumulator(mg1, flow_director="Steepest")
    fa1.run_one_step()

    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x**2 + mg2.node_y**2, at="node")
    fa2 = LossyFlowAccumulator(mg2, flow_director=FlowDirectorSteepest)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x**2 + mg3.node_y**2, at="node")
    fd = FlowDirectorSteepest(mg3)
    fa3 = LossyFlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key], mg1.at_node[key])

        assert_array_equal(mg1.at_node[key], mg2.at_node[key])

        assert_array_equal(mg2.at_node[key], mg3.at_node[key])


def test_director_adding_methods_are_equivalent_D8():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg0.add_field("topographic__elevation", mg0.node_x**2 + mg0.node_y**2, at="node")
    fa0 = LossyFlowAccumulator(mg0, flow_director="D8")
    fa0.run_one_step()

    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x**2 + mg1.node_y**2, at="node")
    fa1 = LossyFlowAccumulator(mg1, flow_director="FlowDirectorD8")
    fa1.run_one_step()

    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x**2 + mg2.node_y**2, at="node")
    fa2 = LossyFlowAccumulator(mg2, flow_director=FlowDirectorD8)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x**2 + mg3.node_y**2, at="node")
    fd = FlowDirectorD8(mg3)
    fa3 = LossyFlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key], mg1.at_node[key])

        assert_array_equal(mg1.at_node[key], mg2.at_node[key])

        assert_array_equal(mg2.at_node[key], mg3.at_node[key])


def test_director_adding_methods_are_equivalent_Dinf():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg0.add_field("topographic__elevation", mg0.node_x**2 + mg0.node_y**2, at="node")
    fa0 = LossyFlowAccumulator(mg0, flow_director="DINF")
    fa0.run_one_step()

    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x**2 + mg1.node_y**2, at="node")
    fa1 = LossyFlowAccumulator(mg1, flow_director="FlowDirectorDINF")
    fa1.run_one_step()

    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x**2 + mg2.node_y**2, at="node")
    fa2 = LossyFlowAccumulator(mg2, flow_director=FlowDirectorDINF)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x**2 + mg3.node_y**2, at="node")
    fd = FlowDirectorDINF(mg3)
    fa3 = LossyFlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key], mg1.at_node[key])

        assert_array_equal(mg1.at_node[key], mg2.at_node[key])

        assert_array_equal(mg2.at_node[key], mg3.at_node[key])


def test_director_adding_methods_are_equivalent_MFD():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg0.add_field("topographic__elevation", mg0.node_x**2 + mg0.node_y**2, at="node")
    fa0 = LossyFlowAccumulator(mg0, flow_director="MFD")
    fa0.run_one_step()

    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x**2 + mg1.node_y**2, at="node")
    fa1 = LossyFlowAccumulator(mg1, flow_director="FlowDirectorMFD")
    fa1.run_one_step()

    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x**2 + mg2.node_y**2, at="node")
    fa2 = LossyFlowAccumulator(mg2, flow_director=FlowDirectorMFD)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg3.add_field("topographic__elevation", mg3.node_x**2 + mg3.node_y**2, at="node")
    fd = FlowDirectorMFD(mg3)
    fa3 = LossyFlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key], mg1.at_node[key])

        assert_array_equal(mg1.at_node[key], mg2.at_node[key])

        assert_array_equal(mg2.at_node[key], mg3.at_node[key])


def test_passing_a_bad_component():
    """Check that a random component can't be a director."""
    from landlab.components import ChiFinder

    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, "topographic__elevation", flow_director=ChiFinder)


def test_error_for_to_many_with_depression():
    """Check that an error is thrown when to_many methods started DF."""

    mg0 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg0.add_field("topographic__elevation", mg0.node_x**2 + mg0.node_y**2, at="node")

    mg1 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg1.add_field("topographic__elevation", mg1.node_x**2 + mg1.node_y**2, at="node")

    with pytest.raises(NotImplementedError):
        LossyFlowAccumulator(
            mg0, flow_director="MFD", depression_finder="DepressionFinderAndRouter"
        )
    with pytest.raises(NotImplementedError):
        LossyFlowAccumulator(
            mg0, flow_director="DINF", depression_finder="DepressionFinderAndRouter"
        )

    fa0 = LossyFlowAccumulator(mg0, flow_director="MFD")
    fa0.run_one_step()
    with pytest.raises(NotImplementedError):
        DepressionFinderAndRouter(mg0)

    fa1 = LossyFlowAccumulator(mg1, flow_director="DINF")
    fa1.run_one_step()
    with pytest.raises(NotImplementedError):
        DepressionFinderAndRouter(mg1)


def test_fields():
    """Check to make sure the right fields have been created.

    Check that the sizes are also correct.
    """
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = LossyFlowAccumulator(mg)
    fa.run_one_step()

    assert sorted(mg.at_node) == [
        "drainage_area",
        "flow__data_structure_delta",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__sink_flag",
        "flow__upstream_node_order",
        "surface_water__discharge",
        "surface_water__discharge_loss",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]

    mg2 = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg2.add_field("topographic__elevation", mg2.node_x + mg2.node_y, at="node")
    fa2 = LossyFlowAccumulator(mg2, flow_director="MFD")
    fa2.run_one_step()
    assert sorted(mg2.at_node.keys()) == [
        "drainage_area",
        "flow__data_structure_delta",
        "flow__link_to_receiver_node",
        "flow__receiver_node",
        "flow__receiver_proportions",
        "flow__sink_flag",
        "flow__upstream_node_order",
        "surface_water__discharge",
        "surface_water__discharge_loss",
        "topographic__elevation",
        "topographic__steepest_slope",
        "water__unit_flux_in",
    ]


@pytest.mark.parametrize("fd", ["Steepest", "D8", "MFD", "DINF"])
def test_accumulated_area_closes(fd):
    """Check that accumulated area is area of core nodes."""
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = LossyFlowAccumulator(mg)
    fa.run_one_step()

    drainage_area = mg.at_node["drainage_area"]
    drained_area = np.sum(drainage_area[mg.boundary_nodes])
    core_area = np.sum(mg.cell_area_at_node[mg.core_nodes])
    assert drained_area == core_area


def test_specifying_routing_method_wrong():
    """Test specifying incorrect method for routing compatability with
    DepressionFinderAndRouter.
    """
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    with pytest.raises(ValueError):
        LossyFlowAccumulator(
            mg, flow_director="D4", depression_finder="DepressionFinderAndRouter"
        )

    df = DepressionFinderAndRouter(mg)
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director="D4", depression_finder=df)


def test_field_name_array_float_case1():
    """Topography as field, runoff rate as float"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    mg.add_field("topographic__elevation", topographic__elevation, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, "topographic__elevation", runoff_rate=10.0)
    assert_array_equal(
        mg.at_node["water__unit_flux_in"], 10.0 * np.ones(mg.size("node"))
    )

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            10.0,
            50.0,
            0.0,
            0.0,
            10.0,
            50.0,
            0.0,
            0.0,
            10.0,
            40.0,
            0.0,
            0.0,
            10.0,
            20.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_field_name_array_float_case2():
    """Topography as field, runoff rate as field name"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    runoff_rate = [
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0,
        4.0,
        5.0,
        5.0,
        5.0,
        5.0,
    ]

    mg.add_field("topographic__elevation", topographic__elevation, at="node")
    mg.add_field("runoff_rate", runoff_rate, at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, "topographic__elevation", runoff_rate="runoff_rate")

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            2.0,
            16.0,
            0.0,  # KRB double checked these numbers by hand 5/15/18 - OK
            0.0,
            2.0,
            16.0,
            0.0,
            0.0,
            3.0,
            14.0,
            0.0,
            0.0,
            4.0,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_field_name_array_float_case3():
    """Topography as field, runoff rate as array"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    runoff_rate = [
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0,
        4.0,
        5.0,
        5.0,
        5.0,
        5.0,
    ]

    mg.add_field("topographic__elevation", topographic__elevation, at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, "topographic__elevation", runoff_rate=runoff_rate)

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            2.0,
            16.0,
            0.0,  # KRB double checked these numbers by hand 5/15/18 - OK
            0.0,
            2.0,
            16.0,
            0.0,
            0.0,
            3.0,
            14.0,
            0.0,
            0.0,
            4.0,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_field_name_array_float_case4():
    """Topography as array, runoff rate as float"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    mg.add_field("topographic__elevation", topographic__elevation, at="node")
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, topographic__elevation, runoff_rate=10.0)
    assert_array_equal(
        mg.at_node["water__unit_flux_in"], 10.0 * np.ones(mg.size("node"))
    )

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            10.0,
            50.0,
            0.0,
            0.0,
            10.0,
            50.0,
            0.0,
            0.0,
            10.0,
            40.0,
            0.0,
            0.0,
            10.0,
            20.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_field_name_array_float_case5():
    """Topography as array, runoff rate as field name"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    runoff_rate = [
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0,
        4.0,
        5.0,
        5.0,
        5.0,
        5.0,
    ]

    mg.add_field("topographic__elevation", topographic__elevation, at="node")
    mg.add_field("runoff_rate", runoff_rate, at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, "topographic__elevation", runoff_rate="runoff_rate")

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            2.0,
            16.0,
            0.0,  # KRB double checked these numbers by hand 5/15/18 - OK
            0.0,
            2.0,
            16.0,
            0.0,
            0.0,
            3.0,
            14.0,
            0.0,
            0.0,
            4.0,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_field_name_array_float_case6():
    """Topography as array, runoff rate as array"""
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
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
            20.0,
            0.0,
            0.0,
            32.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    runoff_rate = [
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0,
        4.0,
        5.0,
        5.0,
        5.0,
        5.0,
    ]

    mg.add_field("topographic__elevation", topographic__elevation, at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = LossyFlowAccumulator(mg, topographic__elevation, runoff_rate=runoff_rate)

    fa.run_one_step()
    reciever = np.array(
        [0, 1, 2, 3, 4, 1, 2, 7, 8, 10, 6, 11, 12, 14, 10, 15, 16, 17, 18, 19]
    )

    da = np.array(
        [
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            5.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.0,
            0.0,
            1.0,
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    q = np.array(
        [
            0.0,
            2.0,
            16.0,
            0.0,  # KRB double checked these numbers by hand 5/15/18 - OK
            0.0,
            2.0,
            16.0,
            0.0,
            0.0,
            3.0,
            14.0,
            0.0,
            0.0,
            4.0,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(mg.at_node["flow__receiver_node"], reciever)
    assert_array_equal(mg.at_node["drainage_area"], da)
    assert_array_equal(mg.at_node["surface_water__discharge"], q)


def test_flow_accumulator_properties():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = LossyFlowAccumulator(mg)
    fa.run_one_step()

    node_drainage_area = np.array(
        [
            0.0,
            3.0,
            3.0,
            3.0,
            0.0,
            0.0,
            3.0,
            3.0,
            3.0,
            0.0,
            0.0,
            2.0,
            2.0,
            2.0,
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

    node_order_upstream = np.array(
        [
            0,
            1,
            6,
            11,
            16,
            2,
            7,
            12,
            17,
            3,
            8,
            13,
            18,
            4,
            5,
            9,
            10,
            14,
            15,
            19,
            20,
            21,
            22,
            23,
            24,
        ]
    )

    assert_array_equal(fa.node_order_upstream, node_order_upstream)
    assert_array_equal(fa.node_water_discharge, node_drainage_area)
    assert_array_equal(fa.node_drainage_area, node_drainage_area)


def test_bad_director_name():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director="spam")


def test_bad_director_instance():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    ld = LinearDiffuser(mg, linear_diffusivity=1.0)
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director=ld)


def test_instantiated_director_with_kwargs():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fd = FlowDirectorSteepest(mg)
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director=fd, partition_method="eggs")


def test_depression_finder_as_bad_string():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director="D8", depression_finder="spam")


def test_depression_finder_as_string():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    LossyFlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )


def test_depression_finder_as_instance():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    df = DepressionFinderAndRouter(mg)
    LossyFlowAccumulator(mg, flow_director="D8", depression_finder=df)


def test_depression_finder_bad_instance():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    ld = LinearDiffuser(mg, linear_diffusivity=1.0)
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director="D8", depression_finder=ld)


def test_instantiated_depression_finder_with_kwargs():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    df = DepressionFinderAndRouter(mg)
    with pytest.raises(ValueError):
        LossyFlowAccumulator(
            mg, flow_director="D8", depression_finder=df, routing="eggs"
        )


def test_depression_finder_bad_uninstantiated_component():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    with pytest.raises(ValueError):
        LossyFlowAccumulator(mg, flow_director="D8", depression_finder=LinearDiffuser)


def test_hex_mfd():
    mg = HexModelGrid((5, 3))
    mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")
    fa = LossyFlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()
