import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import NormalFault


def test_dx_equals_zero():
    """Test a vertical fault trace."""
    grid = RasterModelGrid((6, 6), xy_spacing=10)

    grid.add_zeros("topographic__elevation", at="node")

    param_dict = {
        "faulted_surface": "topographic__elevation",
        "fault_dip_angle": 90.0,
        "fault_throw_rate_through_time": {"time": [0, 9, 10], "rate": [0, 0, 0.05]},
        "fault_trace": {"y1": 0, "x1": 30, "y2": 30, "x2": 30},
        "include_boundaries": True,
    }

    nf = NormalFault(grid, **param_dict)

    out = np.array(
        [
            [True, True, True, False, False, False],
            [True, True, True, False, False, False],
            [True, True, True, False, False, False],
            [True, True, True, False, False, False],
            [True, True, True, False, False, False],
            [True, True, True, False, False, False],
        ],
        dtype=bool,
    )

    assert_array_equal(nf.faulted_nodes.reshape(grid.shape), out)


def test_anti_aximuth_greq_2pi():
    """Test anti azimuth over 2*pi."""
    grid = RasterModelGrid((6, 6), xy_spacing=10)

    grid.add_zeros("topographic__elevation", at="node")

    param_dict = {
        "faulted_surface": "topographic__elevation",
        "fault_dip_angle": 90.0,
        "fault_throw_rate_through_time": {"time": [0, 9, 10], "rate": [0, 0, 0.05]},
        "fault_trace": {"y1": 30.0, "x1": 30.0, "y2": 20.0, "x2": 0.0},
        "include_boundaries": True,
    }

    nf = NormalFault(grid, **param_dict)

    assert nf._fault_anti_azimuth > 2.0 * np.pi

    out = np.array(
        [
            [True, True, True, True, True, True],
            [True, True, True, True, True, True],
            [True, True, True, True, True, True],
            [False, False, False, False, True, True],
            [False, False, False, False, False, False],
            [False, False, False, False, False, False],
        ],
        dtype=bool,
    )

    assert_array_equal(nf.faulted_nodes.reshape(grid.shape), out)


def test_non_raster():
    """Test a hex model grid."""
    grid = HexModelGrid((7, 3), spacing=10, xy_of_lower_left=(-15.0, 0.0))

    grid.add_zeros("topographic__elevation", at="node")

    param_dict = {
        "faulted_surface": "topographic__elevation",
        "fault_dip_angle": 90.0,
        "fault_throw_rate_through_time": {"time": [0, 9, 10], "rate": [0, 0, 0.05]},
        "fault_trace": {"y1": 30.0, "x1": 30.0, "y2": 20.0, "x2": 0.0},
        "include_boundaries": True,
    }

    nf = NormalFault(grid, **param_dict)

    # plotting, to test this. it works!
    # import matplotlib.pyplot as plt
    # plt.figure()
    # imshow_grid(grid, nf.faulted_nodes, color_for_background='y')
    # plt.plot(grid.x_of_node, grid.y_of_node, 'c.')
    # plt.plot([param_dict['fault_trace']['x1'], param_dict['fault_trace']['x2']],
    #         [param_dict['fault_trace']['y1'], param_dict['fault_trace']['y2']], 'r')
    # plt.show()

    out = np.array(
        [
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            False,
            True,
            True,
            True,
            True,
            False,
            False,
            False,
            False,
            True,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
        ],
        dtype=bool,
    )

    assert_array_equal(nf.faulted_nodes, out)


def test_dip_geq_90():
    """Test dip angles of >90 degrees."""
    grid = RasterModelGrid((6, 6), xy_spacing=10)

    grid.add_zeros("topographic__elevation", at="node")

    with pytest.raises(ValueError):
        NormalFault(grid, fault_dip_angle=90.001)


def test_uplifting_multiple_fields():
    """Test uplifting multiple fields with NormalFault."""
    grid = RasterModelGrid((6, 6), xy_spacing=10)

    grid.add_zeros("topographic__elevation", at="node")

    zbr = grid.add_zeros("bedrock__elevation", at="node")

    zbr -= 1.0

    param_dict = {
        "faulted_surface": ["topographic__elevation", "bedrock__elevation"],
        "fault_dip_angle": 90.0,
        "fault_throw_rate_through_time": {"time": [0, 9, 10], "rate": [1, 1, 1]},
        "fault_trace": {"y1": 30.0, "x1": 30.0, "y2": 20.0, "x2": 0.0},
        "include_boundaries": True,
    }

    nf = NormalFault(grid, **param_dict)

    nf.run_one_step(dt=10)

    elev = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            10.0,
            10.0,
            10.0,
            0.0,
            0.0,
            10.0,
            10.0,
            10.0,
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    bedrock = np.array(
        [
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            9.0,
            9.0,
            9.0,
            9.0,
            -1.0,
            -1.0,
            9.0,
            9.0,
            9.0,
            9.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            9.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
        ]
    )

    assert_array_equal(grid.at_node["topographic__elevation"], elev)
    assert_array_equal(grid.at_node["bedrock__elevation"], bedrock)


def test_uplifting_a_not_yet_created_field():
    """Test uplifting a field that does not exist with  NormalFault."""
    grid = RasterModelGrid((6, 6), xy_spacing=10)

    grid.add_zeros("topographic__elevation", at="node")

    zbr = grid.add_zeros("bedrock__elevation", at="node")

    zbr -= 1.0

    param_dict = {
        "faulted_surface": [
            "topographic__elevation",
            "bedrock__elevation",
            "spam",
            "eggs",
        ],
        "fault_dip_angle": 90.0,
        "fault_throw_rate_through_time": {"time": [0, 9, 10], "rate": [1, 1, 1]},
        "fault_trace": {"y1": 30.0, "x1": 30.0, "y2": 20.0, "x2": 0.0},
        "include_boundaries": True,
    }

    assert "spam" not in grid.at_node
    assert "eggs" not in grid.at_node

    # instantiating NormalFault will not create spam or eggs
    nf = NormalFault(grid, **param_dict)

    assert "spam" not in grid.at_node
    assert "eggs" not in grid.at_node

    assert "spam" in nf._not_yet_instantiated
    assert "eggs" in nf._not_yet_instantiated

    # running NormalFault will not create spam or eggs
    nf.run_one_step(dt=10)

    assert "spam" not in grid.at_node
    assert "eggs" not in grid.at_node

    assert "spam" in nf._not_yet_instantiated
    assert "eggs" in nf._not_yet_instantiated

    # running NormalFault after adding spam and eggs will result in NormalFault
    # modifying these fields.

    grid.add_zeros("eggs", at="node")

    grid.add_zeros("spam", at="node")

    nf.run_one_step(dt=10)

    assert "spam" in grid.at_node
    assert "eggs" in grid.at_node

    vals = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            10.0,
            10.0,
            10.0,
            0.0,
            0.0,
            10.0,
            10.0,
            10.0,
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_equal(grid.at_node["eggs"], vals)
    assert_array_equal(grid.at_node["spam"], vals)
