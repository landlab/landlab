import numpy as np
import pytest
from numpy.testing import assert_array_equal
from six import StringIO

from landlab import (
    CLOSED_BOUNDARY,
    HexModelGrid,
    NetworkModelGrid,
    RadialModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
)
from landlab.grid.hex import from_dict as hex_from_dict
from landlab.grid.raster import from_dict as raster_from_dict


def test_raster_old_from_dict_deprecated():
    params = {"NUM_COLS": 10, "NUM_ROWS": 5, "GRID_SPACING": 4}
    with pytest.deprecated_call():
        mg = raster_from_dict(params)
    assert mg.shape == (5, 10)
    assert mg.dx == 4
    assert mg.dy == 4


def test_hex_old_from_dict_deprecated():
    params = {"NUM_COLS": 10, "NUM_ROWS": 5, "GRID_SPACING": 4}
    with pytest.deprecated_call():
        mg = hex_from_dict(params)
    assert mg.number_of_nodes == 54


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
        "axis_name:\n"
        "    - 'spam'\n"
        "    - 'eggs'\n"
        "axis_units:\n"
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
    assert np.all(mg.status_at_node[mg.boundary_nodes] == CLOSED_BOUNDARY)
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
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    mg = RasterModelGrid.from_dict(params)

    # assert things.
    assert mg.shape == (10, 20)
    assert mg.dx == 25
    assert mg.dy == 45
    assert (mg.x_of_node.min(), mg.y_of_node.min()) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes] == CLOSED_BOUNDARY)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_hex_from_dict():
    params = {
        "base_num_rows": 5,
        "base_num_cols": 4,
        "dx": 2.0,
        "xy_of_lower_left": (35, 55),
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
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
        "num_shells": 5,
        "dr": 2.0,
        "xy_of_center": (35, 55),
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
        "xy_of_reference": (12345, 678910),
    }

    mg = RadialModelGrid.from_dict(params)

    # assert things.
    assert mg.number_of_nodes == 95
    assert mg.xy_of_center == (35, 55)
    assert [35, 55] in mg.xy_of_node
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_network_from_dict():
    params = {
        "yx_of_node": [(0, 1, 2, 2), (0, 0, -1, 1)],
        "links": ((1, 0), (2, 1), (3, 1)),
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
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
    file_strn = (
        "yx_of_node:\n"
        "    - [0, 1, 2, 2]\n"
        "    - [0, 0, -1, 1]\n"
        "links:\n"
        "    - [1, 0]\n"
        "    - [2, 1]\n"
        "    - [3, 1]\n"
        "xy_of_reference:\n"
        "    - 12345\n"
        "    - 678910\n"
        "axis_name:\n"
        "    - 'spam'\n"
        "    - 'eggs'\n"
        "axis_units:\n"
        "    - 'smoot'\n"
        "    - 'parsec'"
    )
    file_like = StringIO(file_strn)
    mg = NetworkModelGrid.from_file(file_like)

    assert_array_equal(mg.x_of_node, np.array([0.0, 0.0, -1.0, 1.0]))
    assert_array_equal(mg.y_of_node, np.array([0.0, 1.0, 2.0, 2.0]))
    assert_array_equal(mg.nodes_at_link, np.array([[0, 1], [2, 1], [1, 3]]))
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")
    assert mg.xy_of_reference == (12345, 678910)


def test_voronoi_from_dict():
    x = [0, 0.1, 0.2, 0.3, 1, 1.1, 1.2, 1.3, 2, 2.1, 2.2, 2.3]
    y = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    params = {
        "x": x,
        "y": y,
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
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
    assert_array_equal(mg.node_x, true_x)
    assert_array_equal(mg.node_y, true_y)
    assert_array_equal(mg.adjacent_nodes_at_node, true_nodes_at_node)
