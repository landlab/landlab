import numpy as np
import pytest

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorSteepest
from landlab.utils.distance_to_divide import calculate_distance_to_divide


def test_no_flow_receivers():
    """Test that correct error is raised when no flow recievers are
    on the grid."""
    mg = RasterModelGrid((30, 70))
    with pytest.raises(FieldError):
        calculate_distance_to_divide(mg)


def test_no_upstream_array():
    """Test that correct error is raised when no flow__upstream_node_order."""
    mg = RasterModelGrid((30, 70))
    mg.add_ones("topographic__elevation", at="node")
    mg.add_ones("drainage_area", at="node")
    fd = FlowDirectorSteepest(mg)
    fd.run_one_step()
    with pytest.raises(FieldError):
        calculate_distance_to_divide(mg)


def test_drainage_area():
    """Test that correct error is raised when no flow__upstream_node_order."""
    mg = RasterModelGrid((30, 70))
    mg.add_ones("topographic__elevation", at="node")
    mg.add_ones("flow__upstream_node_order", at="node")
    fd = FlowDirectorSteepest(mg)
    fd.run_one_step()
    with pytest.raises(FieldError):
        calculate_distance_to_divide(mg)


@pytest.mark.parametrize("flow_dir", ["D8", "D4", "MFD"])
def test_simple_case_same(flow_dir):
    mg = RasterModelGrid((20, 5))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.y_of_node
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fa = FlowAccumulator(mg, flow_director=flow_dir)
    fa.run_one_step()

    dist1 = calculate_distance_to_divide(mg, longest_path=True)
    dist2 = calculate_distance_to_divide(mg, longest_path=False)

    np.testing.assert_array_equal(dist1, dist2)

    correct = np.array(
        [
            [0.0, 18.0, 18.0, 18.0, 0.0],
            [0.0, 17.0, 17.0, 17.0, 0.0],
            [0.0, 16.0, 16.0, 16.0, 0.0],
            [0.0, 15.0, 15.0, 15.0, 0.0],
            [0.0, 14.0, 14.0, 14.0, 0.0],
            [0.0, 13.0, 13.0, 13.0, 0.0],
            [0.0, 12.0, 12.0, 12.0, 0.0],
            [0.0, 11.0, 11.0, 11.0, 0.0],
            [0.0, 10.0, 10.0, 10.0, 0.0],
            [0.0, 9.0, 9.0, 9.0, 0.0],
            [0.0, 8.0, 8.0, 8.0, 0.0],
            [0.0, 7.0, 7.0, 7.0, 0.0],
            [0.0, 6.0, 6.0, 6.0, 0.0],
            [0.0, 5.0, 5.0, 5.0, 0.0],
            [0.0, 4.0, 4.0, 4.0, 0.0],
            [0.0, 3.0, 3.0, 3.0, 0.0],
            [0.0, 2.0, 2.0, 2.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    np.testing.assert_array_equal(dist1.reshape(mg.shape), correct)


def test_complex_case():
    mg = RasterModelGrid((20, 5))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.y_of_node
    middle = mg.x_of_node == 2
    z[middle] *= 0.1

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fa = FlowAccumulator(mg, flow_director="D4")
    fa.run_one_step()

    long_path = calculate_distance_to_divide(mg, longest_path=True)
    short_path = calculate_distance_to_divide(mg, longest_path=False)

    long_correct = np.array(
        [
            [0.0, 1.0, 19.0, 1.0, 0.0],
            [0.0, 0.0, 18.0, 0.0, 0.0],
            [0.0, 0.0, 17.0, 0.0, 0.0],
            [0.0, 0.0, 16.0, 0.0, 0.0],
            [0.0, 0.0, 15.0, 0.0, 0.0],
            [0.0, 0.0, 14.0, 0.0, 0.0],
            [0.0, 0.0, 13.0, 0.0, 0.0],
            [0.0, 0.0, 12.0, 0.0, 0.0],
            [0.0, 0.0, 11.0, 0.0, 0.0],
            [0.0, 0.0, 10.0, 0.0, 0.0],
            [0.0, 0.0, 9.0, 0.0, 0.0],
            [0.0, 0.0, 8.0, 0.0, 0.0],
            [0.0, 0.0, 7.0, 0.0, 0.0],
            [0.0, 0.0, 6.0, 0.0, 0.0],
            [0.0, 0.0, 5.0, 0.0, 0.0],
            [0.0, 0.0, 4.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    short_correct = np.array(
        [
            [0.0, 1.0, 3.0, 1.0, 0.0],
            [0.0, 0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    np.testing.assert_array_equal(short_path.reshape(mg.shape), short_correct)
    np.testing.assert_array_equal(long_path.reshape(mg.shape), long_correct)
