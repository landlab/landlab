from io import StringIO

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import (
    HexModelGrid,
    NetworkModelGrid,
    RadialModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
)


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
        [
            37.0,
            39.0,
            41.0,
            43.0,
            36.0,
            38.0,
            40.0,
            42.0,
            44.0,
            35.0,
            37.0,
            39.0,
            41.0,
            43.0,
            45.0,
            36.0,
            38.0,
            40.0,
            42.0,
            44.0,
            37.0,
            39.0,
            41.0,
            43.0,
        ]
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
