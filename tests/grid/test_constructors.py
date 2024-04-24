from io import StringIO

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import FramedVoronoiGrid
from landlab import HexModelGrid
from landlab import NetworkModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid


def test_raster_from_file():
    file_strn = (
        "shape:\n"
        "    - 10\n"
        "    - 20\n"
        "xy_spacing:\n"
        "    - 25\n"
        "    - 45\n"
        "bc:\n"
        "    right: 'closed'\n"
        "    top: 'closed'\n"
        "    left: 'closed'\n"
        "    bottom: 'closed'\n"
        "xy_of_reference:\n"
        "    - 12345\n"
        "    - 678910\n"
        "xy_of_lower_left:\n"
        "    - 35\n"
        "    - 55\n"
        "xy_axis_name:\n"
        "    - 'spam'\n"
        "    - 'eggs'\n"
        "xy_axis_units:\n"
        "    - 'smoot'\n"
        "    - 'parsec'"
    )
    file_like = StringIO(file_strn)
    mg = RasterModelGrid.from_file(file_like)

    # assert things.
    assert mg.shape == (10, 20)
    assert mg.dx == 25
    assert mg.dy == 45
    assert (mg.x_of_node.min(), mg.y_of_node.min()) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes] == mg.BC_NODE_IS_CLOSED)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_raster_from_dict():
    params = {
        "shape": (10, 20),
        "xy_spacing": (25, 45),
        "bc": {
            "right": "closed",
            "top": "closed",
            "left": "closed",
            "bottom": "closed",
        },
        "xy_of_lower_left": (35, 55),
        "xy_axis_name": ("spam", "eggs"),
        "xy_axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    mg = RasterModelGrid.from_dict(params)

    # assert things.
    assert mg.shape == (10, 20)
    assert mg.dx == 25
    assert mg.dy == 45
    assert (mg.x_of_node.min(), mg.y_of_node.min()) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes] == mg.BC_NODE_IS_CLOSED)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_hex_from_dict():
    params = {
        "shape": (5, 4),
        "spacing": 2.0,
        "xy_of_lower_left": (35, 55),
        "xy_axis_name": ("spam", "eggs"),
        "xy_axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    mg = HexModelGrid.from_dict(params)

    # assert things.
    true_x_node = np.array(
        [37.0, 39.0, 41.0, 43.0, 36.0, 38.0, 40.0, 42.0]
        + [44.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 36.0]
        + [38.0, 40.0, 42.0, 44.0, 37.0, 39.0, 41.0, 43.0]
    )
    assert_array_equal(true_x_node, mg.x_of_node)
    assert (mg.x_of_node.min(), mg.y_of_node.min()) == (35, 55)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_radial_from_dict():
    params = {
        "n_rings": 5,
        "nodes_in_first_ring": 4,
        "xy_of_center": (35, 55),
        "xy_axis_name": ("spam", "eggs"),
        "xy_axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    grid_a = RadialModelGrid.from_dict(params)
    grid_b = RadialModelGrid(**params)

    assert grid_a.number_of_nodes == grid_b.number_of_nodes
    assert grid_a.xy_of_center == grid_b.xy_of_center == (35, 55)
    assert grid_a.xy_of_node == pytest.approx(grid_b.xy_of_node)
    assert grid_a.axis_units == ("smoot", "parsec")
    assert grid_a.axis_name == ("spam", "eggs")
    assert grid_a.xy_of_reference == (12345, 678910)


def test_network_from_dict():
    params = {
        "yx_of_node": [(0, 1, 2, 2), (0, 0, -1, 1)],
        "links": ((1, 0), (2, 1), (3, 1)),
        "xy_axis_name": ("spam", "eggs"),
        "xy_axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }
    mg = NetworkModelGrid.from_dict(params)
    assert_array_equal(mg.x_of_node, np.array([0.0, 0.0, -1.0, 1.0]))
    assert_array_equal(mg.y_of_node, np.array([0.0, 1.0, 2.0, 2.0]))
    assert_array_equal(mg.nodes_at_link, np.array([[0, 1], [2, 1], [1, 3]]))
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_network_from_file():
    file_like = StringIO(
        """
        yx_of_node:
            - [0, 1, 2, 2]
            - [0, 0, -1, 1]
        links:
            - [1, 0]
            - [2, 1]
            - [3, 1]
        xy_of_reference:
            - 12345
            - 678910
        xy_axis_name:
            - 'spam'
            - 'eggs'
        xy_axis_units:
            - 'smoot'
            - 'parsec'
        """
    )
    mg = NetworkModelGrid.from_file(file_like)

    assert_array_equal(mg.x_of_node, np.array([0.0, 0.0, -1.0, 1.0]))
    assert_array_equal(mg.y_of_node, np.array([0.0, 1.0, 2.0, 2.0]))
    assert_array_equal(mg.nodes_at_link, np.array([[0, 1], [2, 1], [1, 3]]))
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_framed_voronoi_edge_nodes():
    grid = FramedVoronoiGrid(
        (5, 6), xy_spacing=(10, 15), xy_of_lower_left=(1, 2), xy_min_spacing=7.5
    )

    expected_x = (
        np.array(
            [
                [0, 10, 20, 30, 40, 50],
                [0, 99, 99, 99, 99, 50],
                [0, 99, 99, 99, 99, 50],
                [0, 99, 99, 99, 99, 50],
                [0, 10, 20, 30, 40, 50],
            ]
        )
        + 1.0
    )
    expected_y = (
        np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [15 - 7.5 / 2 - 0.001, 99, 99, 99, 99, 15 + 7.5 / 2 + 0.001],
                [30 - 7.5 / 2 - 0.001, 99, 99, 99, 99, 30 + 7.5 / 2 + 0.001],
                [45 - 7.5 / 2 - 0.001, 99, 99, 99, 99, 45 + 7.5 / 2 + 0.001],
                [60, 60, 60, 60, 60, 60],
            ]
        )
        + 2.0
    )
    actual_x = grid.x_of_node.reshape(grid.shape)
    actual_y = grid.y_of_node.reshape(grid.shape)

    assert_array_almost_equal(expected_x[0, :], actual_x[0, :])
    assert_array_almost_equal(expected_x[-1, :], actual_x[-1, :])
    assert_array_almost_equal(expected_x[:, 0], actual_x[:, 0])
    assert_array_almost_equal(expected_x[:, -1], actual_x[:, -1])

    assert_array_almost_equal(expected_y[0, :], actual_y[0, :])
    assert_array_almost_equal(expected_y[-1, :], actual_y[-1, :])
    assert_array_almost_equal(expected_y[:, 0], actual_y[:, 0])
    assert_array_almost_equal(expected_y[:, -1], actual_y[:, -1])


@pytest.mark.parametrize("xy_min_spacing", (1, 0, 3, (3, 4), (1, 10)))
def test_framed_voronoi_min_spacing(xy_min_spacing):
    grid = FramedVoronoiGrid(
        (100, 200), xy_spacing=(3, 10), xy_min_spacing=xy_min_spacing
    )

    assert np.all(grid.length_of_link >= np.min(xy_min_spacing))


@pytest.mark.parametrize("xy_min_spacing", (10, (3, 10), (4, 11)))
def test_framed_voronoi_bad_min_spacing(xy_min_spacing):
    with pytest.raises(ValueError):
        FramedVoronoiGrid((10, 20), xy_spacing=(3, 10), xy_min_spacing=xy_min_spacing)


def test_framed_voronoi_from_dict():
    params = {
        "shape": (5, 6),
        "xy_spacing": (10.0, 15.0),
        "xy_of_lower_left": (1.0, 2.0),
        "xy_min_spacing": (5.0, 7.5),
        "seed": 200,
        "xy_of_reference": (0.5, 3.0),
        "xy_axis_name": ("a", "b"),
        "xy_axis_units": "l",
    }

    mg = FramedVoronoiGrid.from_dict(params)
    assert mg.shape == (5, 6)
    assert mg.xy_spacing == (10.0, 15.0)
    assert mg._xy_of_lower_left == (1.0, 2.0)
    assert mg._xy_min_spacing == (5.0, 7.5)
    assert mg._seed == 200
    assert mg.xy_of_reference == (0.5, 3.0)
    assert mg.axis_name == ("a", "b")
    assert mg.axis_units == ("l", "l")

    true_x = np.zeros(30)
    true_y = np.zeros(30)
    true_x[0:7] = np.array([1.0, 11.0, 21.0, 31.0, 41.0, 51.0, 1.0])
    true_x[11:13] = np.array([51.0, 1.0])
    true_x[17:19] = np.array([51.0, 1.0])
    true_x[23:30] = np.array([51.0, 1.0, 11.0, 21.0, 31.0, 41.0, 51.0])
    true_y[0:7] = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 13.249])
    true_y[11:13] = np.array([20.751, 28.249])
    true_y[17:19] = np.array([35.751, 43.249])
    true_y[23:30] = np.array([50.751, 62.0, 62.0, 62.0, 62.0, 62.0, 62.0])

    for i in [0, 1, 2, 3, 4, 5, 6, 11, 12, 17, 18, 23, 24, 25, 26, 27, 28, 29]:
        assert_array_almost_equal(mg.x_of_node[i], true_x[i])
        assert_array_almost_equal(mg.y_of_node[i], true_y[i])


def test_voronoi_from_dict():
    x = [0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3]
    y = [0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0]
    params = {
        "x": x,
        "y": y,
        "xy_axis_name": ("spam", "eggs"),
        "xy_axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    mg = VoronoiDelaunayGrid.from_dict(params)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)

    true_x = np.array([0.0, 1.0, 2.0, 0.1, 1.1, 2.1, 0.2, 1.2, 2.2, 0.3, 1.3, 2.3])
    true_y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0])
    true_nodes_at_node = np.array(
        [
            [1, 3, -1, -1, -1, -1],
            [2, 4, 3, 0, -1, -1],
            [5, 4, 1, -1, -1, -1],
            [4, 6, 0, 1, -1, -1],
            [5, 7, 6, 3, 1, 2],
            [8, 7, 4, 2, -1, -1],
            [7, 9, 3, 4, -1, -1],
            [8, 10, 9, 6, 4, 5],
            [11, 10, 7, 5, -1, -1],
            [10, 6, 7, -1, -1, -1],
            [11, 9, 7, 8, -1, -1],
            [10, 8, -1, -1, -1, -1],
        ]
    )
    assert_array_equal(mg.x_of_node, true_x)
    assert_array_equal(mg.y_of_node, true_y)
    assert_array_equal(mg.adjacent_nodes_at_node, true_nodes_at_node)
