import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.field.errors import FieldError
from landlab.grid.raster_mappers import map_max_of_link_nodes_to_link


def test_max_at_link():
    grid = RasterModelGrid((3, 4))

    z = [
        [0, 1, 2, 3],
        [3, 4, 5, 6],
        [1, 2, 3, 4],
    ]
    expected = [
        *[1, 2, 3],
        *[3, 4, 5, 6],
        *[4, 5, 6],
        *[3, 4, 5, 6],
        *[2, 3, 4],
    ]

    actual = map_max_of_link_nodes_to_link(grid, z)
    assert_array_equal(actual, expected)

    grid.at_node["foo"] = z
    actual = map_max_of_link_nodes_to_link(grid, "foo")
    assert_array_equal(actual, expected)


def test_max_at_link_with_out():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("foo", at="node")
    out = np.empty(grid.number_of_links)

    actual = map_max_of_link_nodes_to_link(grid, "foo", out=out)
    assert actual is out
    assert_array_equal(actual, 1.0)


def test_max_at_link_missing_field():
    grid = RasterModelGrid((3, 4))
    with pytest.raises(FieldError, match="foo"):
        map_max_of_link_nodes_to_link(grid, "foo")
